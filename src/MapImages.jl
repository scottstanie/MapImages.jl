__precompile__()


module MapImages

export MapImage

using Parameters
import Sario
import Sario: DemRsc  # Heavily used here
import Base: size, similar, parent, getindex, setindex!



# Note: modelling off OffsetArrays
# https://github.com/JuliaArrays/OffsetArrays.jl/blob/master/src/OffsetArrays.jl
# mutable struct MapImage{T, N, AA<:AbstractArray} <: AbstractArray{T, N}
@with_kw mutable struct MapImage{T, N} <: AbstractArray{T, N}
    image::AbstractArray{T, N}
    demrsc::DemRsc
    # filename::Union{Nothing, AbstractString}  # Do I care to store this??
end
MapImage(A::AbstractArray{T,N}, demrsc::DemRsc) where {T,N} =
    MapImage{T,N}(A, demrsc)


function _get_dem_h5(filename)
    if "dem_rsc" in names(filename)
        demrsc = Sario.load_dem_from_h5(filename)
    else
        demrscfile = Sario.find_rsc_file(filename)
        demrsc = Sario.load(demrscfile)
    end
    return demrsc
end

# TODO: add permutes to pass through to load
function MapImage(filename::AbstractString, dem_or_dset::AbstractString)
    if Sario.get_file_ext(filename) == ".h5"
        demrsc = _get_dem_h5(filename)
        return MapImage(Sario.load(filename, dset_name=dem_or_dset), demrsc)
    else
        demrscfile = Sario.find_rsc_file(filename)
        return MapImage(Sario.load(filename), Sario.load(demrscfile))
    end
end

function MapImage(filename::AbstractString)
    if Sario.get_file_ext(filename) == ".h5"
        demrsc = _get_dem_h5(filename)
        return MapImage(Sario.load(filename), demrsc)
    else
        demrscfile = Sario.find_rsc_file(filename)
        return MapImage(Sario.load(filename), Sario.load(demrscfile))
    end
end
#
# TODO: other loading strategies? or do we want to restrict

Base.eachindex(::IndexCartesian, A::MapImage) = CartesianIndices(axes(A))

# Base.parent(A::MapImage) = A.image  # TODO: do we need?
Base.size(A::MapImage) = size(A.image)
Base.size(A::MapImage, d) = size(A.image, d)

function Base.similar(A::MapImage, ::Type{T}, dims::Dims) where T
    B = similar(A.image, T, dims)
    return MapImage(B, copy(A.demrsc))
end

Base.parent(A::MapImage) = A.image

# Catchall: just getindex on image, ignore demrsc
Base.getindex(A::MapImage{T,N}, I::Vararg{Any,N}) where {T,N} = A.image[I...]

# For single number getting:
# Base.getindex(A::MapImage{T,N}, I::Vararg{Int,N}) where {T,N} = A.image[I...]

ColonOrRange = Union{Colon, <:AbstractRange}
step(::Colon) = 1  # For DemRsc adjustment purposes
start(x::ColonOrRange) = typeof(x) <: AbstractRange ? x.start : 1
length(x::ColonOrRange) = typeof(x) <: AbstractRange ? length(x) : nothing

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

# TODO: tests
# TODO: 3D
# TODO: floats for lat/lon (do i want this?)

# function Base.setindex!(A::MapImage{T,N}, val, I::Vararg{Int,N}) where {T,N}
function Base.setindex!(A::MapImage{T,N}, val, I::Vararg{Any,N}) where {T,N}
    # setindex!(A, X, inds...)
    # A[inds...] = X
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
