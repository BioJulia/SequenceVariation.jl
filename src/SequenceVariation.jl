module SequenceVariation

# TODO: Add functionality to move a Variation to a new reference
# needs to be done with heavy checks to make sure the alignment of the two 
# is not full of gaps in the area where the Variation is
#=
Substitution should work just if the ref nuc don't match a deletion
Insertions/deletions should work if there are only matches (+/- 3 nucs), or at the ends

Perhaps first make an array of oldref -> newref, with some positions marked nothing?
=#
# TODO: Do we need to prevent calling indels at the ends of alignments?
# No, we need a filtering step that trims the alignment before calling Variation

#=
Extract alignment
Trim ends of gaps
Call Variations

=#

# TODO: Add inversions?

using BioSequences
using BioAlignments
import BioAlignments: OP_START, OP_SEQ_MATCH, OP_SEQ_MISMATCH, OP_INSERT, OP_DELETE
import BioSymbols: BioSymbol

abstract type Edit end

struct Substitution{T <: BioSymbol} <: Edit
    symbol::T
end

struct Deletion <: Edit
    len::Int
end

struct Insertion{S <: BioSequence} <: Edit
    seq::S
end

Base.:(==)(x::Insertion{A}, y::Insertion{A}) where A = x.seq == y.seq

struct Diff{E <: Edit}
    pos::Int
    edit::E
end

struct Variant{S <: BioSequence, E <: Edit}
    ref::S
    diffs::Vector{Diff{E}}
end

struct Variation{S <: BioSequence, E <: Edit}
    ref::S
    diff::Diff{E}
end

###
function check(v::Variation{S, <:Substitution{T}}) where {S, T}
    T == eltype(S) || throw(TypeError(:check, "", eltype(S), T))
    checkbounds(v.ref, v.diff.pos)
end

check(v::Variation{S, Deletion}) where S = checkbounds(v.ref, v.diff.pos:(v.diff.pos+v.diff.edit.len)-1)

function check(v::Variation{S, <:Insertion}) where S
    length(v.diff.edit.seq) > 0 || throw(ArgumentError("Insertions cannot be length 0"))
    # We can have insertions at the very end, after the reference sequence
    v.diff.pos == lastindex(v.ref) + 1 && return nothing
    checkbounds(v.ref, v.diff.pos)
end

Base.show(io::IO, x::Diff{<:Substitution}) = print(io, x.pos, x.edit.symbol)
Base.show(io::IO, x::Diff{Deletion}) = print(io, 'Δ', x.pos, '-', x.pos + x.edit.len - 1)
Base.show(io::IO, x::Diff{<:Insertion}) = print(io, x.pos, x.edit.seq)

function Base.show(io::IO, x::Variant{S, <:Substitution}) where S
    print(io, x.ref[x.diff.pos], x.diff.pos, x.diff.edit.symbol)
end

Base.show(io::IO, x::Variant{S, Deletion}) where S = show(io, x.diff)
Base.show(io::IO, x::Variant{S, <:Insertion} where S) = show(io, x.diff)

Base.:(==)(x::T, y::T) where {T <: Variant} = (x.ref === y.ref) & (x.diff == y.diff)

#################

#=
function variations(ref::S, refaln::S, seqaln::S) where {S <: BioSequence}
    aln = AlignedSequence(seqaln, refaln)
    return variations(ref, refaln, seqaln, aln.aln.anchors)
end

function variations(ref::S, refaln::S, seqaln::S, anchors::Vector{AlignmentAnchor}) where {S <: BioSequence{A}} where A
    result = SeqVar{A}[]
    isempty(anchors) && return result
    firstop = first(anchors)
    if (firstop.op !== OP_START) || (firstop.refpos != 0 | firstop.seqpos != 0)
        throw(ArgumentError("Alignment must begin with OP_START at positions (0, 0)"))
    end
    # refpos and seqpos refer to the ungapped sequences. seqaln refer to the gapped one
    refpos, seqpos, seqalnpos, seqalnend = 1, 1, 1, 1
    for anchor in @view anchors[2:end]
        seqend, refend = anchor.seqpos, anchor.refpos
        seqalnend = seqalnpos + max(seqend - seqpos, refend - refpos)

        # For mismatches, we add one Variation per mismatch
        if anchor.op == OP_SEQ_MISMATCH
            p = seqalnpos
            for pos in refpos:refend
                push!(result, SeqVar(ref, pos, Substitution{eltype(A)}(seqaln[p])))
                p += 1
            end

        # Deletions are trivial
        elseif anchor.op == OP_DELETE
            len = refend - refpos + 1
            push!(result, SeqVar(ref, refpos, Deletion(len)))

        # Insertions are a little more tricky
        elseif anchor.op == OP_INSERT
            push!(result, SeqVar(ref, refend + 1, Insertion{A}(seqaln[seqalnpos:seqalnend])))
        
        elseif anchor.op == OP_SEQ_MATCH
            nothing
        else
            throw(ArgumentError("Cannot create Variation from operation $(anchor.op)"))
        end
        refpos, seqpos, seqalnpos = refend + 1, seqend + 1, seqalnend + 1
    end
    return result
end

function posmap(gap, oldref::BioSequence, newref::BioSequence)
    if length(oldref) != length(newref)
        throw(ArgumentError("Sequences must be equal-length and aligned"))
    end
    result = Vector{Union{Int, Nothing}}(undef, length(oldref))
    oldpos, newpos = 0, 0
    for (o, n) in zip(oldref, newref)
        oldpos += o != gap
        newpos += n != gap
        o != gap && (result[oldpos] = n == gap ? nothing : newpos)
    end
    return resize!(result, oldpos)
end

"""
    posmap(oldref::S, newref::S)

Given two aligned sequences, creates a vector `v` such that `v[o] = n`, where `o` is
a position in the old (ungapped) sequence and `n` is a position in the new ungapped sequence.
If `o` maps to a gap in the new sequence, `v[o] === nothing`.
"""
posmap(oldref::S, newref::S) where {A, S <: BioSequence{A}} = posmap(gap(eltype(A)), oldref, newref)

function rereference(var::SeqVar{A, S}, posmap, ref::LongSequence{A}) where {A, S <: Substitution}
    newpos = posmap[var.pos]
    newpos === nothing && throw(ArgumentError("Position $(var.pos) maps to a gap in reference"))
    return SeqVar{A, S}(ref, newpos::Int, var.var)
end

function checkflanks(posmap, varpos::Int, flank::Int)
    for i in max(1, varpos - flank) : min(length(posmap), varpos + flank)
        posmap[i] === nothing && throw(ArgumentError("Position $i maps to a deletion"))
    end
end

function rereference(var::SeqVar{A, Deletion}, posmap, ref::LongSequence{A}, flank::Int=5) where A
    checkbounds(posmap, var.pos:var.pos + var.var.len - 1)
    checkflanks(posmap, var.pos, flank)
    return SeqVar{A, Deletion}(ref, posmap[var.pos], var.var)
end

function rereference(var::SeqVar{A, Insertion{A}}, posmap, ref::LongSequence{A}, flank::Int=5) where A
    checkflanks(posmap, var.pos, flank)
    return SeqVar{A, Insertion{A}}(ref, posmap[var.pos], var.var)
end


# Assumes same reference, and that they are not mutally exclusive
# (e.g no substitution in deleted areas)
#=
function reconstruct(v::Vector{<:SeqVar{A}}) where A
    isempty(v) && throw(ArgumentError("Need at least one Variation to reconstruct sequence"))
    srt = sort(v, by=x -> x.pos)
    len = length(v[1].ref)
    for i in srt
        if i isa SeqVar{A, <:Deletion}
            len -= i.var.len
        elseif i isa SeqVar{A, <:Insertion}
            len += length(i.var.seq)
        end
    end


    ref = copy(v[1].ref)
    oldpos, newpos = 1, 1
    for i in srt
=#
=#


end # module
