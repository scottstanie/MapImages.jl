__precompile__()


module MapImages

export MapImage


import Parameters: @with_kw, @unpack
import Sario
import Sario: DemRsc
import Base: size, similar, step, parent, getindex, setindex!


# TODO: tests

# Note: modelling off OffsetArrays
# https://github.com/JuliaArrays/OffsetArrays.jl/blob/master/src/OffsetArrays.jl
# mutable struct MapImage{T, N, AA<:AbstractArray} <: AbstractArray{T, N}
@with_kw mutable struct MapImage{T, N} <: AbstractArray{T, N}
    image::AbstractArray{T, N}
    demrsc::DemRsc
    # filename::Union{Nothing, AbstractString}  # Do I care to store this??
end
MapImage(A::AbstractArray{T,N}, demrsc::DemRsc) where {T,N} = MapImage{T,N}(A, demrsc)


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
function MapImage(filename::AbstractString, dem_or_dset::AbstractString)
    demrsc = _get_dem(filename)
    return MapImage(Sario.load(filename, dset_name=dem_or_dset), demrsc)
end

function MapImage(filename::AbstractString)
    demrsc = _get_dem(filename)
    return MapImage(Sario.load(filename), demrsc)
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
Base.getindex(A::MapImage{T,N}, I::Vararg{Any,N}) where {T,N} = A.image[I...]


ColonOrRange = Union{Colon, <:AbstractRange}
# For DemRsc adjustment purposes
Base.step(x::Colon) = 1
start(x::ColonOrRange) = typeof(x) <: AbstractRange ? x.start : 1
length(x::ColonOrRange) = typeof(x) <: AbstractRange ? Base.length(x) : nothing

# 2D array subarray: handles ranges and colons
function Base.getindex(A::MapImage{T,2}, I::Vararg{Union{Colon, <:AbstractRange}, 2}) where {T,N}
    subimg = A.image[I...]
    
    # DemRsc adjusting:
    rowrange, colrange = I
    newdemrsc = _new_dem_data(A.demrsc, rowrange, colrange)

    return MapImage(subimg, newdemrsc)
end

function _new_dem_data(d::DemRsc, rowrange::ColonOrRange, colrange::ColonOrRange)
    new_rows = isnothing(length(rowrange)) ? d.rows : length(rowrange)
    new_cols = isnothing(length(colrange)) ? d.cols : length(colrange)
    row_step = step(rowrange)
    col_step = step(colrange)
    row_start = start(rowrange)
    col_start = start(colrange)
    return crop_rsc_data(d, new_rows, new_cols,
                         row_start, col_start,
                         row_step, col_step)
end

"""Adjusts the old demrsc of a MapImage after slicing

Takes the 'rows' and 'cols' from demrsc
and adjusts for the smaller size with a new dict

The row/col_start are the indices where we have sliced:
    e.g. img[2:3:11, :] would have 
    row_start = 2, row_step = 3, new_rows = 4
"""
function crop_rsc_data(demrsc::DemRsc, 
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
                  x_step=new_x_step, y_step=new_y_step)
end

# TODO: 3D
# TODO: floats for lat/lon (do i want this?)

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


# Extra functions manipulating MapImages
include("./latlon.jl")

end # module
