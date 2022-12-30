"""
    Insertion{S <: BioSequence}

Represents the insertion of a `S` into a sequence. The location of the insertion
is stored outside the struct.
"""
struct Insertion{S<:BioSequence}
    seq::S

    function Insertion{S}(x::S) where {S<:BioSequence}
        isempty(x) && error("Insertion must be at least 1 symbol")
        return new(x)
    end
end
Insertion(s::BioSequence) = Insertion{typeof(s)}(s)

Base.length(x::Insertion) = length(x.seq)
Base.:(==)(x::Insertion, y::Insertion) = x.seq == y.seq
Base.hash(x::Insertion, h::UInt) = hash(Insertion, hash(x.seq, h))

function _refbases(i::Insertion, reference::S, pos::UInt) where {S<:BioSequence}
    return S([reference[pos]])
end

function _altbases(i::Insertion, reference::S, pos::UInt) where {S<:BioSequence}
    if pos == 1
        return S([i.seq..., reference[pos]])
    else
        return S([reference[pos], i.seq...])
    end
end
