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
@with_kw mutable struct MapImage{T, N, AA<:AbstractArray} <: AbstractArray{T, N}
    image::AA
    demrsc::DemRsc
    # filename::Union{Nothing, AbstractString}  # Do I care to store this??
end
MapImage(A::AbstractArray{T,N}, demrsc::DemRsc) where {T,N} =
    MapImage{T,N,typeof(A)}(A, demrsc)

function MapImage(filename::AbstractString, dem_or_dset::AbstractString)
    if Sario.get_file_ext(filename) == ".h5"
        demrsc = Sario.load_dem_from_h5(filename)
        return MapImage(Sario.load(filename, dset_name=dem_or_dset), demrsc)
    else
        demrscfile = Sario.find_rsc_file(filename)
        return MapImage(Sario.load(filename), Sario.load(demrscfile))
    end
end

function MapImage(filename::AbstractString)
    if Sario.get_file_ext(filename) == ".h5"
        demrsc = Sario.load_dem_from_h5(filename)
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

function Base.getindex(A::MapImage{T,N}, I::Vararg{Int,N}) where {T,N}
    ret = A.image[I...]
    # TODO: add the demrsc adjustments
end

function Base.setindex!(A::MapImage{T,N}, val, I::Vararg{Int,N}) where {T,N}
    # setindex!(A, X, inds...)
    # A[inds...] = X
    A.image[I...] = val
    val
end


#TODO:
#https://docs.julialang.org/en/v1/base/io-network/index.html#Base.show-Tuple{Any,Any,Any}
# MIME"image/png"
# show(io, ::MIME"image/png", x::MapImage) = .


# Extra functions manipulating MapImages
include("./latlon.jl")

end # module
