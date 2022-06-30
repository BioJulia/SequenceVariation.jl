module SequenceVariation

"""
Needs to be able to:
* Given a sequence and a reference, create a `Variant` that unambiguously represents
the sequence

* Given a `Variant` and a new reference, translate the variant to the new reference.

* Given a mutation and a reference and a sequence, determine if the sequence has that
mutation

TODO now:
* Play around with some NGS results rel. to picked reference.
    * Is it easy to construct ref and variants? I.e. is API nice?
    * Is it nice and easy to check if a mut is present?
    *

* Implement "reference switching".
* Add tests
"""

using BioAlignments: BioAlignments, PairwiseAlignment
using BioGenerics: BioGenerics, leftposition
using BioSequences: BioSequences, BioSequence, NucleotideSeq, AminoAcidSeq, LongSequence, isgap
using BioSymbols: BioSymbol

const BA = BioAlignments
const BS = BioSequences

#=
import Automa
import Automa.RegExp: @re_str
=#

struct Unsafe end
struct Inapplicable end

#=
const CONTEXT = Automa.CodeGenContext(
    vars=Automa.Variables(:p, :p_end, :p_eof, :ts, :te, :cs, :data, :mem, :byte),
    generator=:goto,
    checkbounds=false
)


const BYTES = Union{String, SubString{String}, Vector{UInt8}}
=#

"""
    Substitution

Represents the presence of a `T` at a given position. The position is stored
outside this struct.
"""
struct Substitution{T <: BioSymbol}
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
        new(len)
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
struct Insertion{S <: BioSequence}
    seq::S

    function Insertion{S}(x::S) where {S <: BioSequence}
        isempty(x) && error("Insertion must be at least 1 symbol")
        new(x)
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
struct Edit{S <: BioSequence, T <: BioSymbol}
    x::Union{Substitution{T}, Deletion, Insertion{S}}
    pos::UInt
end
Base.:(==)(e1::Edit, e2::Edit) = e1.pos == e2.pos && e1.x == e2.x
Base.hash(x::Edit, h::UInt) = hash(Edit, hash((x.x, x.pos), h))

function Base.parse(::Type{T}, s::AbstractString) where {T <: Edit{Se, Sy}}  where {Se, Sy}
    parse(T, String(s))
end

function Base.parse(::Type{<:Edit{Se, Sy}}, s::Union{String, SubString{String}}) where {Se, Sy}
    # Either "Δ1-2", "11T" or "G16C"
    if (m = match(r"^Δ(\d+)-(\d+)$", s); m) !== nothing
        pos = parse(UInt, m[1])
        stop = parse(UInt, m[2])
        stop ≥ pos || throw(ArgumentError("Non-positive deletion length: \"" * s * "\""))
        Edit{Se, Sy}(Deletion(stop - pos + 1), pos)
    elseif (m = match(r"^(\d+)([A-Za-z]+)$", s); m) !== nothing
        pos = parse(UInt, m[1])
        seq = Se(m[2])
        Edit{Se, Sy}(Insertion(seq), pos)
    elseif (m = match(r"^[A-Za-z](\d+)([A-Za-z])$", s); m) !== nothing
        pos = parse(UInt, m[1])
        sym = Sy(first(m[2]))
        Edit{Se, Sy}(Substitution(sym), pos)
    else
        throw(ArgumentError("Failed to parse edit \"" * s * '"'))
    end
end

mutation(e::Edit) = e.x
BioGenerics.leftposition(e::Edit) = e.pos

#=
@noinline throw_parse_error(T, p::Integer) = error("Failed to parse $T at byte $p")

# Parse substitution
let
    machine = let
        biosymbol = re"[A-Za-z]"
        number = re"[0-9]+"

        biosymbol.actions[:enter] = [:enter]
        number.actions[:all] = [:digit]

        Automa.compile(biosymbol * number * biosymbol)
    end
    actions = Dict(
        :enter => quote
            symbol = T(Char(byte))
            if firstsymbol === nothing
                firstsymbol = symbol
            else
                lastsymbol = symbol
            end
        end,
        :digit => :(num = UInt(10)*num + (byte - 0x30) % UInt),
    )
    @eval begin
        function Base.parse(::Type{Edit{S, T}}, data::BYTES) where {S, T}
            $(Automa.generate_init_code(CONTEXT, machine))
            p_eof = p_end = sizeof(data)
            firstsymbol = lastsymbol = nothing
            num = UInt(0)
            $(Automa.generate_exec_code(CONTEXT, machine, actions))
            iszero(cs) || throw_parse_error(Edit{S, T}, p)
            if firstsymbol == lastsymbol
                error("First symbol and last symbol are identical")
            end
            return Edit{S, T}(Substitution{T}(lastsymbol), num)
        end
    end
end
=#

# Edits are applied sequentially from first to last pos.
# The vector must always be sorted by pos.
struct Variant{S <: BioSequence, T <: BioSymbol}
    ref::S
    edits::Vector{Edit{S, T}}

    Variant{S, T}(ref::S, edits::Vector{Edit{S, T}}, ::Unsafe) where {S, T} = new(ref, edits)
end

function Variant{S,T}(ref::S, edits::Vector{Edit{S, T}}) where {S<:BioSequence, T<:BioSymbol}
    sort!(edits, by=x -> x.pos)
    result = Variant{S, T}(ref, edits, Unsafe())
    is_valid(result) || error("TODO") # report what kind of error message?
    return result
end

function Variant(ref::S, edits::Vector{Edit{S, T}}) where {S<:BioSequence, T<:BioSymbol}
    Variant{S, T}(ref, edits)
end

function Base.show(io::IO, x::Variant)
    n = length(x.edits)
    print(io, summary(x), " with $n edit$(n > 1 ? "s" : ""):")
    for i in x.edits
        v = Variation(x.ref, i)
        print(io, "\n  ")
        show(io, v)
    end
end

# Validate:
# A sequence is invalid if any of its operations are out of bounds, or the same position
# is affected by multiple edits.
function is_valid(v::Variant)
    isempty(v.ref) && return false
    valid_positions = 1:length(v.ref)
    last_was_insert = false
    for edit in v.edits
        pos = edit.pos
        op = edit.x
        # For substitutions we simply do not allow another modification of the same base
        if op isa Substitution
            pos in valid_positions || return false
            valid_positions = first(valid_positions) + 1 : last(valid_positions)
            last_was_insert = false
        # Insertions affect 0 reference bases, so it does not modify the valid positions
        # for next op. However, we cannot have two insertions at the same position, because
        # then the order of them is ambiguous
        elseif op isa Insertion
            pos in (first(valid_positions)-1+last_was_insert:last(valid_positions)) || return false
            last_was_insert = true
        # Deletions obviously invalidate the reference bases that are deleted.
        elseif op isa Deletion
            len = length(op)
            pos in (first(valid_positions):last(valid_positions)-len+1) || return false
            valid_positions = first(valid_positions) + len : last(valid_positions)
            last_was_insert = false
        end
    end
    return true
end

function Variant(aln::PairwiseAlignment{T, T}) where {T <: LongSequence{<:Union{BS.AminoAcidAlphabet, BS.NucleicAcidAlphabet}}}
    ref = aln.b
    E = eltype(typeof(ref))
    edits = Edit{T, E}[]
    result = Variant(ref, edits)
    refpos = first(aln.a.aln.anchors).refpos
    seqpos = first(aln.a.aln.anchors).seqpos
    markpos = 0
    n_gaps = n_ins = 0
    insertion_buffer = E[]
    for (seqi, refi) in aln
        isgap(refi) || (refpos += 1)
        isgap(seqi) || (seqpos += 1)

        # Check for deletions
        if isgap(seqi)
            iszero(n_gaps) && (markpos = refpos)
            n_gaps += 1
        elseif !iszero(n_gaps)
            push!(edits, Edit{T, E}(Deletion(UInt(n_gaps)), UInt(markpos)))
            n_gaps = 0
        end

        # Check for insertions
        if isgap(refi)
            iszero(n_ins) && (markpos = refpos + 1)
            push!(insertion_buffer, seqi)
            n_ins += 1
        elseif !iszero(n_ins)
            seq = T(insertion_buffer)
            push!(edits, Edit{T, E}(Insertion(seq), UInt(markpos)))
            empty!(insertion_buffer)
            n_ins = 0
        end

        # Substitutions
        if !isgap(refi) && !isgap(seqi) && seqi != refi
            push!(edits, Edit{T, E}(Substitution{E}(seqi), UInt(refpos)))
        end
    end

    # Final indel, if applicable
    if !iszero(n_gaps)
        push!(edits, Edit{T, E}(Deletion(UInt(n_gaps)), UInt(markpos)))
    elseif !iszero(n_ins)
        push!(edits, Edit{T, E}(Insertion(T(insertion_buffer)), UInt(markpos)))
    end

    return result
end

function lendiff(edit::Edit)
    x = edit.x
    x isa Substitution ? 0 : (x isa Deletion ? -length(x) : length(x.x))
end

function reconstruct!(seq::S, x::Variant{S}) where S
    len = length(x.ref) + sum(edit -> lendiff(edit), x.edits)
    resize!(seq, len % UInt)
    refpos = seqpos = 1
    for edit in x.edits
        while refpos < edit.pos
            seq[seqpos] = x.ref[refpos]
            refpos += 1
            seqpos += 1
        end
        editx = edit.x
        if editx isa Substitution
            seq[seqpos] = editx.x
            seqpos += 1
            refpos += 1
        elseif editx isa Deletion
            refpos += editx.len
        elseif editx isa Insertion
            for i in editx.x
                seq[seqpos] = i
                seqpos += 1
            end
        end
    end
    while seqpos ≤ length(seq)
        seq[seqpos] = x.ref[refpos]
        refpos += 1
        seqpos += 1
    end
    seq
end

struct Variation{S <: BioSequence, T <: BioSymbol}
    ref::S
    edit::Edit{S, T}

    function Variation{S, T}(ref::S, e::Edit{S, T}, ::Unsafe) where {S <: BioSequence, T <: BioSymbol}
        new(ref, e)
    end
end

function Variation{S, T}(ref::S, e::Edit{S, T}) where {S <: BioSequence, T <: BioSymbol}
    v = Variation{S, T}(ref, e, Unsafe())
    is_valid(v) ? v : throw(ArgumentError("Invalid variant"))
end

Variation(ref::S, edit::Edit{S, T}) where {S, T} = Variation{S, T}(ref, edit)

function Variation(ref::S, edit::AbstractString) where {S<:BioSequence}
    T = eltype(ref)

    e = parse(Edit{S,T}, edit)
    return Variation{S,T}(ref, e)
end

reference(v::Variation) = v.reference
edit(v::Variation) = v.edit
mutation(v::Variation) = mutation(edit(v))
BioGenerics.leftposition(v::Variation) = leftposition(edit(v))
Base.:(==)(x::Variation, y::Variation) = x.ref == y.ref && x.edit == y.edit
Base.hash(x::Variation, h::UInt) = hash(Variation, hash((x.ref, x.edit), h))

function is_valid(v::Variation)
    isempty(v.ref) && return false
    op = v.edit.x
    pos = v.edit.pos
    if op isa Substitution
        return pos in eachindex(v.ref)
    elseif op isa Insertion
        return pos in 0:lastindex(v.ref)
    elseif op isa Deletion
        return pos in 1:(lastindex(v.ref)-length(op) + 1)
    end
end

function Base.show(io::IO, x::Variation)
    content = x.edit.x
    pos = x.edit.pos
    if content isa Substitution
        print(io, x.ref[pos], pos, content.x)
    elseif content isa Deletion
        print(io, 'Δ', pos, '-', pos + content.len - 1)
    elseif content isa Insertion
        print(io, pos, content.seq)
    else
        print(io, pos, content.x)
    end
end

function Base.in(v::Variation, var::Variant)
    if v.ref != var.ref
        error("References must be equal")
    end
    any(v.edit == edit for edit in var.edits)
end

function translate(var::Variation{S, T}, aln::PairwiseAlignment{S, S}) where {S, T}
    kind = var.edit.x
    pos = var.edit.pos
    seq, ref = aln.seq, aln.b

    # Special case: Insertions may have a pos of 0, which cannot be mapped to
    # the seq using ref2seq
    if iszero(pos)
        (s, r), _ = iterate(aln)
        (isgap(s) | isgap(r)) && return Inapplicable()
        return Variation{S, T}(seq, Edit{S, T}(Insertion(var.edit.x), 0))
    end

    (seqpos, op) = BA.ref2seq(aln, pos)
    if kind isa Substitution
        # If it's a substitution, return nothing if it maps to a deleted
        # position, or substitutes to same base.
        op in (BA.OP_MATCH, BA.OP_SEQ_MATCH, BA.OP_SEQ_MISMATCH) || return nothing
        seq[seqpos] == kind.x && return nothing
        edit = Edit{S, T}(kind, seqpos)
        return Variation{S, T}(seq, edit, Unsafe())
    elseif kind isa Deletion
        # If it's a deletion, return nothing if the deleted part is already missing
        # from the new reference.
        (stop, op2) = BA.ref2seq(aln, pos + length(kind) - 1)
        start = seqpos + op == BA.OP_DELETE
        start < stop && return nothing
        edit = Edit{S, T}(Deletion(stop - start + 1), start)
        return Variation{S, T}(seq, edit, Unsafe())
    else
        # If it maps directly to a symbol, just insert
        if op in (BA.OP_MATCH, BA.OP_SEQ_MATCH, BA.OP_SEQ_MISMATCH)
            # This happens if there is already an insertion at the position
            if pos != lastindex(ref) && first(ref2seq(aln, pos+1)) != seqpos + 1
                return Inapplicable()
            else
                edit = Edit{S, T}(Insertion(var.edit.x), seqpos)
                return Variation{S, T}(seq, edit, Unsafe())
            end
        # Alternatively, it can map to a deletion. In that case, it become really
        # tricky to talk about the "same" insertion.
        else
            return Inapplicable()
        end
    end
end

export Insertion,
    Deletion,
    Substitution,
    Variant,
    Variation,
    reference,
    mutation

end # module
