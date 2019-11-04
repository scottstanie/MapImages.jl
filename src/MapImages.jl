__precompile__(true)



module MapImages

export MapImage

import Parameters: @with_kw, @unpack
import Sario
import Sario: DemRsc
import Base: length, size, similar, step, parent, getindex, setindex!


# Note: modelling off OffsetArrays
# https://github.com/JuliaArrays/OffsetArrays.jl/blob/master/src/OffsetArrays.jl
# mutable struct MapImage{T, N, AA<:AbstractArray} <: AbstractArray{T, N}
@with_kw mutable struct MapImage{T, N} <: AbstractArray{T, N}
    image::AbstractArray{T, N}
    demrsc::DemRsc
    # filename::Union{Nothing, AbstractString}  # Do I care to store this??
end
MapImage(A::AbstractArray{T,N}, demrsc::DemRsc) where {T,N} = MapImage{T,N}(A, demrsc)


# Extra functions manipulating MapImages
include("./latlon.jl")
include("./gps.jl")
include("./geotiff.jl")


# TODO: maybe move this to Sario
function _get_dem(filename)
    if Sario._is_h5(filename) && "dem_rsc" in names(filename)
        demrsc = Sario.load_dem_from_h5(filename)
    else
        demrscfile = Sario.find_rsc_file(filename)
        demrsc = Sario.load(demrscfile)
    end
    return demrsc
end

# TODO: add permutes to pass through to load
function MapImage(filename::AbstractString; dset_name::AbstractString="")
    demrsc = _get_dem(filename)
    return MapImage(Sario.load(filename, dset_name=dset_name), demrsc)
end

Base.eachindex(::IndexCartesian, A::MapImage) = CartesianIndices(axes(A))

Base.parent(A::MapImage) = A.image

# Base.parent(A::MapImage) = A.image  # TODO: do we need?
Base.size(A::MapImage) = size(A.image)
Base.size(A::MapImage, d) = size(A.image, d)

# Taken from https://docs.julialang.org/en/v1/manual/interfaces/ and OffsetArray examples
Base.similar(A::MapImage) = MapImage(similar(A.image), copy(A.demrsc))
Base.similar(A::MapImage, ::Type{S}) where S = MapImage(similar(A.image, S, size(A)), copy(A.demrsc))
Base.similar(A::MapImage, dims::Dims) = MapImage(similar(A.image, eltype(A.image), size(A)),copy(A.demrsc))

function Base.similar(A::MapImage, ::Type{T}, dims::Dims) where T
    B = similar(A.image, T, dims)
    return MapImage(B, copy(A.demrsc))
end


# Allow broadasting, while keeping the output as a MapImage (instead of Array)
# Note: Example taken directly from here:
# https://docs.julialang.org/en/v1/manual/interfaces/#man-interfaces-broadcasting-1
Base.BroadcastStyle(::Type{<:MapImage}) = Broadcast.ArrayStyle{MapImage}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{MapImage}}, ::Type{ElType}) where ElType
    # Scan the inputs for the MapImage:
    A = find_obj(bc)
    # Use the DemRsc field of A to create the output
    MapImage(similar(Array{ElType}, axes(bc)), copy(A.demrsc))
end

"`A = find_obj(As)` returns the first MapImage among the arguments."
find_obj(bc::Base.Broadcast.Broadcasted) = find_obj(bc.args)
find_obj(args::Tuple) = find_obj(find_obj(args[1]), Base.tail(args))
find_obj(x) = x
find_obj(a::MapImage, rest) = a
find_obj(::Any, rest) = find_obj(rest)


# Catchall: just getindex on image, ignore demrsc
function Base.getindex(A::MapImage{T,N}, I::Vararg{Any,N}) where {T,N}
    A.image[I...]
end


ColonOrRange = Union{Colon, <:AbstractRange{Int}}
# For DemRsc adjustment purposes
Base.step(x::Colon) = 1
start(x::ColonOrRange) = typeof(x) <: AbstractRange ? x.start : 1
Base.length(x::ColonOrRange) = typeof(x) <: AbstractRange ? Base.length(x) : nothing

# 2D + 3D array subarray: handles ranges and colons
function Base.getindex(A::MapImage{T,N}, I::Vararg{ColonOrRange, N}) where {T, N}
    subimg = A.image[I...]
    
    # DemRsc adjusting:
    rowrange, colrange = I
    newdemrsc = _new_dem_data(A.demrsc, rowrange, colrange)

    return MapImage(subimg, newdemrsc)
end

# 2D array single floats + colons (for lat/lon)
# First: A[lat, lon]
function Base.getindex(A::MapImage{T,2}, I::Vararg{<:Real, 2}) where {T}
    row = _get_row(A.demrsc, I[1])
    col = _get_col(A.demrsc, I[2])

    subimg = A.image[row, col]
    
    # We assume this is no longer a value image since it's not a range:
    # just return image data
    return subimg
end

# 3D array [lat, lon, :] or [row, col, :]
function Base.getindex(A::MapImage{T,3}, I::Vararg{<:Real, 3}) where {T}
    row = _get_row(A.demrsc, I[1])
    col = _get_col(A.demrsc, I[2])

    # No longer a valid MapImage, just 1d
    return A.image[row, col, I[3]]
end

# 2D array subarray: For float ranges, we want a tuple of floats
# This is because doing something like 30.1:33.4 won't tell you what the
# final number was.
# The StepRangeLen gives the start, len, offset, and step, so you will lose what 
# the end of the range was
# E.G.
#   A[(30.1, 32.2), (-104.1, :)]
function Base.getindex(A::MapImage{T,2}, I::Vararg{Tuple{<:Real, <:Real}, 2}) where {T}
    lats, lons = I

    rowrange, colrange = latlon_to_rowcol(A.demrsc, lats, lons)
    subimg = A.image[rowrange, colrange]

    
    # DemRsc adjusting:
    newdemrsc = _new_dem_data(A.demrsc, rowrange, colrange)

    return MapImage(subimg, newdemrsc)
end


# Making changes to the DemRsc: prep incoming data
function _new_dem_data(d::DemRsc, rowrange::ColonOrRange, colrange::ColonOrRange)
    new_rows = isnothing(length(rowrange)) ? d.rows : length(rowrange)
    new_cols = isnothing(length(colrange)) ? d.cols : length(colrange)
    row_step = step(rowrange)
    col_step = step(colrange)
    row_start = start(rowrange)
    col_start = start(colrange)
    return crop_demrsc(d;
                       new_rows=new_rows, 
                       new_cols=new_cols,
                       row_start=row_start,
                       col_start=col_start,
                       row_step=row_step,
                       col_step=col_step)
end

"""
    crop_demrsc(demrsc::DemRsc;
                new_rows=demrsc.rows, new_cols=demrsc.cols, 
                row_start=1, col_start=1, 
                row_step=1, col_step=1)

Adjusts the old demrsc of a MapImage after slicing.
Takes the 'rows' and 'cols' from demrsc
and adjusts for the smaller size with a new dict

The `row_start`/`col_start` are the indices where we have sliced.
This means:

    img[2:3:11, :] 

would have `row_start` = 2, `row_step` = 3, `new_rows` = 4

"""
function crop_demrsc(demrsc::DemRsc; 
                     new_rows=demrsc.rows, new_cols=demrsc.cols, 
                     row_start=1, col_start=1, 
                     row_step=1, col_step=1)
    @unpack rows, cols, y_step, x_step, y_first, x_first = demrsc

    # Move forward the starting row/col from where it used to be
    new_x_first = x_first + x_step * (col_start - 1)
    new_y_first = y_first + y_step * (row_start - 1)

    # now adjust step sizes (say, if we take only every other row)
    new_x_step = x_step * col_step
    new_y_step = y_step * row_step

    return DemRsc(demrsc, width=new_cols, cols=new_cols,
                  file_length=new_rows, rows=new_rows,
                  x_first=new_x_first, y_first=new_y_first,
                  x_step=new_x_step, y_step=new_y_step)
end



# function Base.setindex!(A::MapImage{T,N}, val, I::Vararg{Int,N}) where {T,N}
function Base.setindex!(A::MapImage{T,N}, val, I::Vararg{Any,N}) where {T,N}
    # setindex!(A, X, inds...)  corresponds to A[inds...] = X 
    A.image[I...] = val
    return val
end

#TODO:
#https://docs.julialang.org/en/v1/base/io-network/index.html#Base.show-Tuple{Any,Any,Any}
# MIME"image/png"
# show(io, ::MIME"image/png", x::MapImage) = .



end # module
