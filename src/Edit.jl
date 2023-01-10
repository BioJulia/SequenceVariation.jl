include("edits/Substitution.jl")
include("edits/Deletion.jl")
include("edits/Insertion.jl")

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

Base.length(e::Edit) = length(_mutation(e))
Base.:(==)(e1::Edit, e2::Edit) = e1.pos == e2.pos && e1.x == e2.x
Base.hash(x::Edit, h::UInt) = hash(Edit, hash((x.x, x.pos), h))
function Base.isless(x::Edit, y::Edit)
    if leftposition(x) == leftposition(y)
        return length(x) < length(y)
    end

    return leftposition(x) < leftposition(y)
end

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

"""
    _mutation(e::Edit)

Returns the underlying [`Substitution`](@ref), [`Insertion`](@ref), or [`Deletion`](@ref) of
`e`.
"""
_mutation(e::Edit) = e.x
BioGenerics.leftposition(e::Edit) = e.pos
function BioGenerics.rightposition(e::Edit)
    if _mutation(e) isa Substitution
        return leftposition(e)
    elseif _mutation(e) isa Insertion
        return leftposition(e) + 1
    elseif _mutation(e) isa Deletion
        return leftposition(e) + length(e) - 1
    else
        error("Unknown mutation type $(typeof(_mutation(e)))")
    end
end

"""
    _lendiff(edit::Edit)

Gets the number of bases that `edit` adds to the reference sequence
"""
function _lendiff(edit::Edit)
    x = _mutation(edit)
    # Each edit type has logic for its length, we just need to know what direction to go
    multiplier = x isa Substitution ? 0 : (x isa Deletion ? -1 : 1)
    return length(x) * multiplier
end
