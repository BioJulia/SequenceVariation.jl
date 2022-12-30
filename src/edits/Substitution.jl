"""
    Substitution

Represents the presence of a `T` at a given position. The position is stored
outside this struct.
"""
struct Substitution{T<:BioSymbol}
    x::T
end

Base.length(::Substitution) = 1
Base.:(==)(x::Substitution, y::Substitution) = x.x == y.x
Base.hash(x::Substitution, h::UInt) = hash(Substitution, hash(x.x, h))

function _refbases(::Substitution, reference::S, pos::UInt) where {S<:BioSequence}
    return S([reference[pos]])
end

function _altbases(s::Substitution, ::S, pos::UInt) where {S<:BioSequence}
    return S([s.x])
end
