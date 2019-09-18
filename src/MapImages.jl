__precompile__()


module MapImages

export MapImage

import Sario

DemRsc = Dict{String, Any}

# Note: modelling off OffsetArrays
mutable struct MapImage{T, N, AA<:AbstractArray} <: AbstractArray{T, N}
    image::AA
    demrsc::DemRsc
    # filename::Union{Nothing, AbstractString}  # Do I care to store this??
end

MapImage(demrscfile::String, filename::String) = MapImage(Sario.load(filename),
                                                          Sario.load(demrscfile),
                                                          filename)
MapImage(demrscfile::String, filename::String) = MapImage(Sario.load(filename),
                                                          Sario.load(demrscfile),
                                                          filename)
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

function Base.getindex(A::MapImage{T,N}, I::Vararg{Int,N}) where {T,N}
    ret = A.image[I...]
    # TODO: add the demrsc adjustments
end

function Base.setindex!(A::OffsetArray{T,N}, val, I::Vararg{Int,N}) where {T,N}
    # setindex!(A, X, inds...)
    # A[inds...] = X
    A.image[I...] = val
    val
end

end # module
