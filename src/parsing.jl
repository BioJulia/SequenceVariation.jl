const DELETION_PATTERN = r"^Δ(\d+)(?:-(\d+))?$"
const INSERTION_PATTERN = r"^(\d+)([A-Z]+)$"
const SUBSTITUTION_PATTERN = r"^([A-Z])(\d+)([A-Z])$"
const MULTI_SUBST_PATTERN = r"^([A-Z](?:/[A-Z])*)(\d+)([A-Z](?:/[A-Z])*)$"

function Base.parse(::Type{<:SeqVar{A}}, ref::LongSequence{A}, st::AbstractString) where A
    m = match(DELETION_PATTERN, st)
    if m !== nothing
        return parse(SeqVar{A, Deletion}, ref, m)
    end
    m = match(SUBSTITUTION_PATTERN, st)
    if m !== nothing
        return parse(SeqVar{A, Substitution{eltype(A)}}, ref, m)
    end
    m = match(INSERTION_PATTERN, st)
    if m !== nothing
        return parse(SeqVar{A, Insertion{A}}, ref, m)
    else
        throw(ArgumentError("Cannot parse $st as Variation"))
    end
end

function Base.parse(::Type{SeqVar{A, Deletion}}, ref::LongSequence{A}, m::RegexMatch) where A
    pos = parse(Int, m[1])
    len = m === nothing ? 1 : parse(Int, m[2]) - pos + 1
    return SeqVar(ref, pos, Deletion(len))
end

function Base.parse(::Type{<:SeqVar{A, Substitution{T}}}, ref::LongSequence{A}, m::RegexMatch) where {A,T}
    if T !== eltype(A)
        throw(ArgumentError("Substitution type must be alphabet eltype"))
    end
    refst, pos, altst = m[1], m[2], m[3]
    pos = parse(Int, pos)
    refT = T(first(refst))
    altT = T(first(altst))
    return SeqVar(ref, pos, Substitution{T}(altT))
end

function Base.parse(::Type{<:SeqVar{A, Insertion{A}}}, ref::LongSequence{A}, m::RegexMatch) where A
    pos = parse(Int, m[1])
    # We accept an insertion immediately after the reference, right?
    checkbounds(ref, pos - (pos == length(ref)))
    seq = LongSequence{A}(m[2])
    return SeqVar(ref, pos, Insertion{A}(seq))
end



totext(m::SeqVar, seqid::IdDict{<:LongSequence, AbstractString}) = string(seqid[m.ref], '\t', m)

function fromtext(st::AbstractString, seqid::Dict{String, LongSequence{A}}) where A
    name, varstr = split(st, '\t', limit=2)
    ref = seqid[name]
    return parse(SeqVar{A}, ref, varstr)
end