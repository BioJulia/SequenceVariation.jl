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

const ALN_MODEL = AffineGapScoreModel(BLOSUM62, gap_open=-12, gap_extend=-2)

@sum_type Edit{S, T}
    Substitution{S, T}(::T)
    Deletion{S, T}(::UInt)
    Insertion{S, T}(::S)
end

"""
    Edit

Abstract type representing a type of nucleotide edit: Deletion, insertion or
substitution.
"""
abstract type Edit end

"""
    Substitution{T <: BioSymbol} <: Edit

Represents the presence of a `T` at a given position. The position is stored
outside this struct.
"""
struct Substitution{T <: BioSymbol} <: Edit
    symbol::T
end

Base.eltype(::Type{<:Substitution{T}}) where T = T

"""
    Deletion <: Edit

Represents the deletion of N symbols. The location of the deletion is stored
outside this struct
"""
struct Deletion <: Edit
    len::UInt
end

"""
    Insertion{S <: BioSequence} <: Edit

Represents the insertion of a `T` into a sequence. The location of the insertion
is stored outside the struct.
"""
struct Insertion{S <: BioSequence} <: Edit
    seq::S
end

Base.:(==)(x::Insertion{A}, y::Insertion{A}) where A = x.seq == y.seq

"""
    Diff{E <: Edit}

Represents an `Edit` og type `E` at a given position.
"""
struct Diff{E <: Edit}
    pos::UInt
    edit::E
end

Diff(pos::Integer, edit::Edit) = Diff{typeof(edit)}(pos, edit)

"""
    Variant{S <: BioSequence, D <: Diff}

Represents a series of diffs of type `D` against a reference of type `S`.
See also `Variation`
"""
struct Variant{S <: BioSequence, D <: Diff}
    ref::S
    diffs::Vector{D}
end

#function Variant(ref::BioSequence, diffs::Vector{D}) where D <: Diff
#    return Variant{typeof(ref), D}(ref, diffs)
#end

function Base.show(io::IO, ::MIME"text/plain", x::Variant{S,D}) where {S,D}
    cols = displaysize(io)[2] - 3
    recur_io = IOContext(io, :SHOWN_SET => x.diffs)
    print(io, summary(x), ":")
    for i in x.diffs
        v = Variation{S, typeof(i)}(x.ref, i)
        str = sprint(show, v, context=recur_io, sizehint=0)
        print(io, "\n  ", Base._truncate_at_width_or_chars(str, cols, "\r\n"))
    end
end

#Variant(s::BioSequence, v::Vector{<:Diff}) = Variant{typeof(s),eltype(v)}(s,v)

"""
    Variation{S <: BioSequence, D <: Diff}

Represent a single diff of type `D` is a sequence of type `S`. See also `Variant`
"""
struct Variation{S <: BioSequence, D <: Diff}
    ref::S
    diff::D
end

#Variation(ref::BioSequence, diff::Diff) = Variation{typeof(ref), typeof(diffs)}(ref, diff)

function check(v::Variation{S, <:Substitution{T}}) where {S, T}
    T == eltype(S) || throw(TypeError(:check, "", eltype(S), T))
    checkbounds(v.ref, v.diff.pos)
end

check(v::Variation{S, Diff{Deletion}}) where S = checkbounds(v.ref, v.diff.pos:(v.diff.pos+v.diff.edit.len)-1)

function check(v::Variation{S, <:Diff{<:Insertion}}) where S
    length(v.diff.edit.seq) > 0 || throw(ArgumentError("Insertions cannot be length 0"))
    # We can have insertions at the very end, after the reference sequence
    v.diff.pos == lastindex(v.ref) + 1 && return nothing
    checkbounds(v.ref, v.diff.pos)
end

Base.show(io::IO, x::Diff{<:Substitution}) = print(io, x.pos, x.edit.symbol)
Base.show(io::IO, x::Diff{Deletion}) = print(io, 'Δ', x.pos, '-', x.pos + x.edit.len - 1)
Base.show(io::IO, x::Diff{<:Insertion}) = print(io, x.pos, x.edit.seq)

function Base.show(io::IO, x::Variation{S, <:Diff{<:Substitution}}) where S
    print(io, x.ref[x.diff.pos], x.diff.pos, x.diff.edit.symbol)
end

Base.show(io::IO, x::Variation{S, Diff{Deletion}}) where S = show(io, x.diff)
Base.show(io::IO, x::Variation{S, <:Diff{<:Insertion}} where S) = show(io, x.diff)

Base.:(==)(x::T, y::T) where {T <: Variation} = (x.ref === y.ref) & (x.diff == y.diff)

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


function variant(ref::S, seq::S) where {S <: BioSequence}
    aln = pairalign(GlobalAlignment(), seq, ref, ALN_MODEL).aln
    diffs = Diff[]
    result = Variant(ref, diffs)
    refpos = 0
    markpos = 0
    n_gaps = n_ins = 0
    for (seqi, refi) in aln
        isgap(refi) || (refpos += 1)

        # Check for deletions
        if isgap(seqi)
            iszero(n_gaps) && (markpos = refi)
            n_gaps += 1
        else
            if !iszero(n_gaps)
                push(diffs, Diff{Deletion}Deletion(n_gaps))
            end
            n_gaps = 0
        end
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
