"""
    Substitution

Represents the presence of a `T` at a given position. The position is stored
outside this struct.
"""
struct Substitution{T<:BioSymbol}
    x::T
end
Base.:(==)(x::Substitution, y::Substitution) = x.x == y.x
Base.hash(x::Substitution, h::UInt) = hash(Substitution, hash(x.x, h))

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

"""
    Edit{S <: BioSequence, T <: BioSymbol}

An edit of either `Substitution{T}`, `Insertion{S}` or `Deletion` at a position.
If deletion: Deletion of length L at ref pos `pos:pos+L-1`
If insertion: Insertion of length L b/w ref pos `pos:pos+1`
"""
struct Edit{S<:BioSequence,T<:BioSymbol}
    x::Union{Substitution{T},Deletion,Insertion{S}}
    pos::UInt
end
Base.:(==)(e1::Edit, e2::Edit) = e1.pos == e2.pos && e1.x == e2.x
Base.hash(x::Edit, h::UInt) = hash(Edit, hash((x.x, x.pos), h))
Base.length(e::Edit) = e isa Substitution ? 1 : length(mutation(e))

function Base.parse(::Type{T}, s::AbstractString) where {T<:Edit{Se,Sy}} where {Se,Sy}
    return parse(T, String(s))
end

function Base.parse(::Type{<:Edit{Se,Sy}}, s::Union{String,SubString{String}}) where {Se,Sy}
    # Either "Δ1-2", "11T" or "G16C"
    if (m = match(r"^Δ(\d+)-(\d+)$", s); m) !== nothing
        pos = parse(UInt, m[1])
        stop = parse(UInt, m[2])
        stop ≥ pos || throw(ArgumentError("Non-positive deletion length: \"" * s * "\""))
        Edit{Se,Sy}(Deletion(stop - pos + 1), pos)
    elseif (m = match(r"^(\d+)([A-Za-z]+)$", s); m) !== nothing
        pos = parse(UInt, m[1])
        seq = Se(m[2])
        Edit{Se,Sy}(Insertion(seq), pos)
    elseif (m = match(r"^[A-Za-z](\d+)([A-Za-z])$", s); m) !== nothing
        pos = parse(UInt, m[1])
        sym = Sy(first(m[2]))
        Edit{Se,Sy}(Substitution(sym), pos)
    else
        throw(ArgumentError("Failed to parse edit \"" * s * '"'))
    end
end

mutation(e::Edit) = e.x
BioGenerics.leftposition(e::Edit) = e.pos
function BioGenerics.rightposition(e::Edit)
    if mutation(e) isa Substitution
        return leftposition(e)
    elseif mutation(e) isa Insertion
        return leftposition(e) + 1
    elseif mutation(e) isa Deletion
        return leftposition(e) + length(e) - 1
    else
        error("Unknown mutation type $(typeof(mutation(e)))")
    end
end

function lendiff(edit::Edit)
    x = edit.x
    return x isa Substitution ? 0 : (x isa Deletion ? -length(x) : length(x.x))
end

function _refbases(s::Substitution, reference::S, pos::UInt) where {S<:BioSequence}
    return S([reference[pos]])
end

function _altbases(s::Substitution, reference::S, pos::UInt) where {S<:BioSequence}
    return S([s.x])
end

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
