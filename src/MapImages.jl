__precompile__()


module MapImages

export MapImage

using Parameters
import Sario
import Sario: DemRsc  # Heavily used here
import Base: size, similar, parent, getindex, setindex!


include("./latlon.jl")

# Note: modelling off OffsetArrays
# https://github.com/JuliaArrays/OffsetArrays.jl/blob/master/src/OffsetArrays.jl
# mutable struct MapImage{T, N, AA<:AbstractArray} <: AbstractArray{T, N}
mutable struct MapImage{T, N, AA<:AbstractArray} <: AbstractArray{T, N}
    image::AA
    demrsc::DemRsc
    # filename::Union{Nothing, AbstractString}  # Do I care to store this??
end
MapImage(A::AbstractArray{T,N}, demrsc::DemRsc) where {T,N} =
    MapImage{T,N,typeof(A)}(A, demrsc)

MapImage(filename::AbstractString, 
         demrscfile::AbstractString) = MapImage(Sario.load(filename),
                                                Sario.load(demrscfile))
function MapImage(filename::AbstractString)
    demrscfile = Sario.find_rsc_file(filename)
    MapImage(Sario.load(filename), Sario.load(demrscfile))
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




end # module
