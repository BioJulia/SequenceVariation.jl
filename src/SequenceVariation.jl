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
import BioSymbols: BioSymbol
using SumTypes

const ALN_MODEL = AffineGapScoreModel(BLOSUM62, gap_open=-12, gap_extend=-2)

@sum_type Edit{S, T} begin
    Substitution{S, T}(::T)
    Deletion{S, T}(::UInt)
    Insertion{S, T}(::S)
end

"""
    Substitution

Represents the presence of a `T` at a given position. The position is stored
outside this struct.
"""
Substitution

Base.eltype(::Type{<:Substitution{S, T}}) where {S, T} = T

"""
    Deletion

Represents the deletion of N symbols. The location of the deletion is stored
outside this struct
"""
Deletion

"""
    Insertion

Represents the insertion of a `T` into a sequence. The location of the insertion
is stored outside the struct.
"""
Insertion

function Insertion(s::BioSequence)
    length(s) == 0 && throw(ArgumentError("Insertions cannot be length 0"))
    Insertion{typeof(s), eltype(s)}(s)
end

Base.:(==)(x::Insertion, y::Insertion) = x.seq == y.seq

"""
    Diff{E <: Edit}

Represents an `Edit` og type `E` at a given position.
"""
struct Diff{S <: BioSequence, T <: BioSymbol}
    pos::UInt
    edit::Edit{S, T}
end

"""
    Variant{S <: BioSequence, T <: BioSymbol}

Represents a series of diffs of type `Diff{S, T}` against a reference of type `S`.
See also `Variation`.
"""
struct Variant{S <: BioSequence, T <: BioSymbol}
    ref::S
    diffs::Vector{Diff{S, T}}
end

function Base.show(io::IO, ::MIME"text/plain", x::Variant{S,T}) where {S,T}
    cols = displaysize(io)[2] - 3
    recur_io = IOContext(io, :SHOWN_SET => x.diffs)
    print(io, summary(x), ":")
    for i in x.diffs
        v = Variation{S, eltype(S)}(x.ref, i)
        str = sprint(show, v, context=recur_io, sizehint=0)
        print(io, "\n  ", Base._truncate_at_width_or_chars(str, cols, "\r\n"))
    end
end

"""
    Variation{S <: BioSequence, T <: BioSymbol}

Represent a single diff against a biosequence. See also `Variant`
"""
struct Variation{S <: BioSequence, T <: BioSymbol}
    ref::S
    diff::Diff{S, T}
end

@case Edit function check((x,)::Substitution{S, T}) where {S, T}
    (ref, pos) -> begin
        eltype(ref) == T || throw(TypeError(:check, "", eltype(ref), T))
        checkbounds(ref, pos)
    end
end

@case Edit function check((seq,)::Insertion{S, T}) where {S, T}
    (ref, pos) -> begin
        # We can have insertions at the very end, after the reference sequence
        pos == lastindex(ref) + 1 && return nothing
        checkbounds(ref, pos)
    end
end

@case Edit function check((len,)::Deletion{S, T}) where {S, T}
    (ref, pos) -> begin
        checkbounds(ref, pos:(pos+len)-1)
    end
end

@assert SumTypes.iscomplete(check, Edit)

function check(x::Variation)
    check(x.diff.edit)(x.ref, x.diff.pos)
end

@case Edit _show((x,)::Substitution) = (io, pos) -> print(io, pos, x)
@case Edit _show((len,)::Deletion) = (io, pos) -> print(io, 'Δ', pos, '-', pos + len - 1)
@case Edit _show((seq,)::Insertion) = (io, pos) -> print(io, pos, seq)

@assert SumTypes.iscomplete(_show, Edit)

Base.show(io::IO, x::Diff) = _show(x.edit)(io, x.pos)

function Base.show(io::IO, x::Variation)
    if x.diff.edit.data isa Substitution
        print(io, x.ref[x.diff.pos])
    end
    show(io, x.diff)
end

Base.:(==)(x::T, y::T) where {T <: Variation} = (x.ref === y.ref) & (x.diff == y.diff)

function variant(seq::LongAminoAcidSeq, ref::LongAminoAcidSeq)
    aln = pairalign(OverlapAlignment(), seq, ref, ALN_MODEL).aln
    diffs = Diff{LongAminoAcidSeq, AminoAcid}[]
    result = Variant(ref, diffs)
    refpos = seqpos = 0
    markpos = 0
    n_gaps = n_ins = 0
    insertion_buffer = AminoAcid[]
    for (seqi, refi) in aln
        isgap(refi) || (refpos += 1)
        isgap(seqi) || (seqpos += 1)

        # Check for deletions
        if isgap(seqi)
            iszero(seqpos) && continue # skip indels at start
            iszero(n_gaps) && (markpos = refpos)
            n_gaps += 1
        else
            if !iszero(n_gaps)
                push!(diffs, Diff(UInt(markpos), Deletion{LongAminoAcidSeq, AminoAcid}(UInt(n_gaps))))
                n_gaps = 0
            end
        end

        # Check for insertions
        if isgap(refi)
            iszero(refpos) && continue # skip indels at start
            iszero(n_ins) && (markpos = refpos + 1)
            push!(insertion_buffer, seqi)
            n_ins += 1
        else
            if !iszero(n_ins)
                seq = LongAminoAcidSeq(insertion_buffer)
                push!(diffs, Diff(UInt(markpos), Insertion(seq)))
                empty!(insertion_buffer)
                n_ins = 0
            end
        end

        # Substitutions
        if !isgap(refi) && !isgap(seqi) && seqi != refi
            push!(diffs, Diff(UInt(refpos), Substitution{LongAminoAcidSeq, AminoAcid}(seqi)))
        end
    end
    # At the end of the loop?
    return result  
end

@case Edit lendiff((x,)::Substitution) = 0
@case Edit lendiff((len,)::Deletion) = -(len % Int)
@case Edit lendiff((seq,)::Insertion) = length(seq) % Int

function reconstruct!(seq::S, x::Variant{S}) where S
    sort!(x.diffs, by=y -> y.pos)
    len = length(x.ref) + sum(diff -> lendiff(diff.edit), x.diffs)
    resize!(seq, len % UInt)
    refpos = seqpos = 1
    for diff in x.diffs
        while refpos < diff.pos
            seq[seqpos] = x.ref[refpos]
            refpos += 1
            seqpos += 1
        end
        if diff.edit.data isa Substitution
            seq[seqpos] = diff.edit.data._1
            seqpos += 1
            refpos += 1
        elseif diff.edit.data isa Deletion
            refpos += diff.edit.data._1
        elseif diff.edit.data isa Insertion
            for i in diff.edit.data._1
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

"""
    reconstruct(x::Variant{S})

Reconstruct the sequence of type `S` that created the variant. It is assumed the
variant is well-formed, e.g. no substitutions in deleted sequences, or
deletions/insertions of the same area multiple times.
"""
reconstruct(x::Variant{S}) where S = reconstruct!(S(), x)

export Substitution, Insertion, Deletion, Diff, Variant, Variation

end # module

#=

function variations(ref::S, refaln::S, seqaln::S) where {S <: BioSequence}
    aln = AlignedSequence(seqaln, refaln)
    return variations(ref, refaln, seqaln, aln.aln.anchors)
end

"Calculate all substitutions between two sequences"
function substitutions(ref::S, seq::S) where {S <: BioSequence}
    if length(ref) != length(seq)
        throw(DimensionMismatch("Sequences must match in length"))
    end
    if !(count(isgap, ref) == count(isgap, seq) == 0)
        throw(ArgumentError("Neither sequence can contain gaps"))
    end
    diffs = Diff{Substitution{eltype(S)}}[]
    @inbounds for i in eachindex(ref)
        refi, seqi = ref[i], seq[i]
        if refi != seqi
            push!(diffs, Diff(i, Substitution(seqi)))
        end
    end
    return Variant(ref, diffs)
end


 
#=
function variations(ref::S, refaln::S, seqaln::S, anchors::Vector{AlignmentAnchor}) where {S <: BioSequence{A}} where A
    result = Variant(dna"TAG", Diff[])
    diffs = result.diffs
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
                push!(diffs, Diff(pos, Substitution{eltype(A)}(seqaln[p])))
                p += 1
            end

        # Deletions are trivial
        elseif anchor.op == OP_DELETE
            len = refend - refpos + 1
            push!(diffs, Diff(refpos, Deletion(len)))

        # Insertions are a little more tricky
        elseif anchor.op == OP_INSERT
            push!(diffs, Diff(refend + 1, Insertion{S}(seqaln[seqalnpos:seqalnend])))
        
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


function rereference(diff::Diff{<:Substitution}, posmap)
    newpos = posmap[diff.pos]
    newpos === nothing && throw(ArgumentError("Position $(diff.pos) maps to a gap in reference"))
    return typeof(diff)(newpos::Int, diff.edit)
end


function checkflanks(posmap, pos::Int, flank::Int)
    for i in max(1, pos - flank) : min(length(posmap), pos + flank)
        posmap[i] === nothing && throw(ArgumentError("Flanking position $i maps to a deletion"))
    end
end

function rereference(diff::Diff{Deletion}, posmap, flank::Int=5)
    checkbounds(posmap, diff.pos:diff.pos + diff.edit.len - 1)
    checkflanks(posmap, diff.pos, flank)
    return typeof(diff)(posmap[diff.pos], diff.edit.len)
end

function rereference(diff::Diff{<:Insertion}, posmap, flank::Int=5)
    checkflanks(posmap, diff.pos, flank)
    return typeof(diff)(posmap[diff.pos], diff.edit)
end


# Assumes same reference, and that they are not mutally exclusive
# (e.g no substitution in deleted areas)
function reconstruct(v::Variant)
    isempty(v.diffs) && return copy(v.ref)
    diffs = sort(v.diffs, by=x -> x.pos)

    # First, get length of result
    len::Int = length(v.ref)
    for diff in diffs
        if diff.edit isa Deletion
            len = len - diff.edit.len
        elseif diff.edit isa Insertion
            len = len + length(diff.edit.seq)
        end
    end

    # Fill in
    dst = typeof(v.ref)(len)
    src, diffi, srci::Int, dsti = v.ref, 1, 1, 1
    diff = diffs[diffi]
    while dsti ≤ len
        while (srci < diff.pos) & (dsti ≤ len)
            dst[dsti] = src[srci]
            srci, dsti = srci + 1, dsti + 1
        end
        diff === nothing && break
        if diff.edit isa Substitution
            dst[dsti] = diff.edit.symbol
            srci, dsti = srci + 1, dsti + 1
        elseif diff.edit isa Insertion
            for i in 1:(length(diff.edit.seq)::Int)
                dst[dsti] = diff.edit.seq[i]
                dsti += 1
            end
        elseif diff.edit isa Deletion
            srci += diff.edit.len
        end
        diffi += 1
        diff = diffi > length(diffs) ? nothing : diffs[diffi]
    end
    return dst
end
=#

export Diff,
    Edit,
    Variation,
    Variant,
    Substitution,
    Insertion,
    Deletion,
    variations,
    substitutions

end # module

=#