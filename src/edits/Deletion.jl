"""
    Deletion

Represents the deletion of N symbols. The location of the deletion is stored
outside this struct
"""
struct Deletion
    len::UInt

    function Deletion(len::UInt)
        iszero(len) && error("Deletion must be at least 1 symbol")
        return new(len)
    end
end
Deletion(x::Integer) = Deletion(convert(UInt, x))

Base.length(x::Deletion) = Int(x.len)
Base.hash(x::Deletion, h::UInt) = hash(Deletion, hash(x.len, h))

function _refbases(d::Deletion, reference::S, pos::UInt) where {S<:BioSequence}
    if pos == 1
        return S(reference[UnitRange{Int}(pos, pos + length(d))])
    else
        return S(reference[UnitRange{Int}(pos - 1, pos + length(d) - 1)])
    end
end

function _altbases(d::Deletion, reference::S, pos::UInt) where {S<:BioSequence}
    if pos == 1
        return S([reference[pos + 1]])
    else
        return S([reference[pos - 1]])
    end
end
