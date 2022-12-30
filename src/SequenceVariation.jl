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

using BioAlignments: BioAlignments, PairwiseAlignment, OP_SOFT_CLIP
using BioGenerics: BioGenerics, leftposition, rightposition
using BioSequences: BioSequences, BioSequence, NucleotideSeq, LongSequence, isgap
using BioSymbols: BioSymbol

const BA = BioAlignments
const BS = BioSequences

struct Unsafe end
struct Inapplicable end

include("Edit.jl")
include("Variant.jl")

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

function Variant(ref::S, vars::Vector{Variation{S,T}}) where {S<:BioSequence, T<:BioSymbol}
    edits = edit.(vars)
    return Variant{S, T}(ref, edits)
end

reference(v::Variation) = v.ref
edit(v::Variation) = v.edit
mutation(v::Variation) = mutation(edit(v))
BioGenerics.leftposition(v::Variation) = leftposition(edit(v))
BioGenerics.rightposition(v::Variation) = rightposition(edit(v))
Base.:(==)(x::Variation, y::Variation) = x.ref == y.ref && x.edit == y.edit
Base.hash(x::Variation, h::UInt) = hash(Variation, hash((x.ref, x.edit), h))

function is_valid(v::Variation)
    isempty(v.ref) && return false
    op = v.edit.x
    pos = v.edit.pos
    if op isa Substitution
        return pos in eachindex(v.ref)
    elseif op isa Insertion
        return pos in 0:lastindex(v.ref)+1
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

function variations(v::Variant{S,T}) where {S,T}
    vs = Vector{Variation{S,T}}(undef, length(edits(v)))
    for (i, e) in enumerate(edits(v))
        vs[i] = Variation{S,T}(reference(v), e)
    end
    return vs
end

function refbases(v::Variation)
    return _refbases(mutation(v), reference(v), leftposition(v))
end

function altbases(v::Variation)
    return _altbases(mutation(v), reference(v), leftposition(v))
end

export Insertion,
    Deletion,
    Substitution,
    Variant,
    Variation,
    reference,
    mutation,
    variations,
    refbases,
    altbases

end # module
