"""
    Haplotype{S<:BioSequence,T<:BioSymbol}

A set of variations within a given sequence that are all found together. Depending on the
field, it might also be referred to as a "genotype" or "strain."

# Constructors

    Haplotype(ref::S, edits::Vector{Edit{S,T}}) where {S<:BioSequence,T<:BioSymbol}
    Haplotype(ref::S, vars::Vector{Variation{S,T}}) where {S<:BioSequence,T<:BioSymbol}
    Haplotype(
        aln::PairwiseAlignment{T,T}
    ) where {T<:LongSequence{<:Union{BS.AminoAcidAlphabet,BS.NucleicAcidAlphabet}}}

When constructing a `Haplotype` from a vector of [`Edit`](@ref)s or [`Variation`](@ref)s,
the edits are applied sequentially from first to last position, therefore the vector must
always be sorted by position. These edits are sorted automatically if constructing from an
alignment.
"""
struct Haplotype{S<:BioSequence,T<:BioSymbol}
    ref::S
    edits::Vector{Edit{S,T}}

    Haplotype{S,T}(ref::S, edits::Vector{Edit{S,T}}, ::Unsafe) where {S,T} = new(ref, edits)
end

function Haplotype{S,T}(
    ref::S, edits::Vector{Edit{S,T}}
) where {S<:BioSequence,T<:BioSymbol}
    sort!(edits; by=x -> x.pos)
    result = Haplotype{S,T}(ref, edits, Unsafe())
    _is_valid(result) || error("TODO") # report what kind of error message?
    return result
end

function Haplotype(ref::S, edits::Vector{Edit{S,T}}) where {S<:BioSequence,T<:BioSymbol}
    return Haplotype{S,T}(ref, edits)
end

function Base.show(io::IO, x::Haplotype)
    n = length(x.edits)
    print(io, summary(x), " with $n edit$(n > 1 ? "s" : ""):")
    for i in x.edits
        v = Variation(x.ref, i)
        print(io, "\n  ")
        show(io, v)
    end
end

"""
    is_valid(h::Haplotype)

Validate `h`. `h` is invalid if any of its operations are out of bounds, or the same
position is affected by multiple edits.
"""
function _is_valid(h::Haplotype)
    isempty(h.ref) && return false
    valid_positions = 1:length(h.ref)
    last_was_insert = false
    for edit in h.edits
        pos = edit.pos
        op = edit.x
        # Sanity check: for this to be a valid variant, it must be comprised of valid
        # variations
        _is_valid(Variation(h.ref, edit)) || return false

        # For substitutions we simply do not allow another modification of the same base
        if op isa Substitution
            pos in valid_positions || return false
            valid_positions = (first(valid_positions) + 1):last(valid_positions)
            last_was_insert = false
            # Insertions affect 0 reference bases, so it does not modify the valid positions
            # for next op. However, we cannot have two insertions at the same position, because
            # then the order of them is ambiguous
        elseif op isa Insertion
            pos in
            ((first(valid_positions) - 1 + last_was_insert):(last(valid_positions) + 1)) ||
                return false
            last_was_insert = true
            # Deletions obviously invalidate the reference bases that are deleted.
        elseif op isa Deletion
            len = length(op)
            pos in (first(valid_positions):(last(valid_positions) - len + 1)) ||
                return false
            valid_positions = (first(valid_positions) + len):last(valid_positions)
            last_was_insert = false
        end
    end
    return true
end

function Haplotype(
    aln::PairwiseAlignment{T,T}
) where {T<:LongSequence{<:Union{BS.AminoAcidAlphabet,BS.NucleicAcidAlphabet}}}
    ref = aln.b
    E = eltype(typeof(ref))
    edits = Edit{T,E}[]
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
            push!(edits, Edit{T,E}(Deletion(UInt(n_gaps)), UInt(markpos)))
            n_gaps = 0
        end

        # Check for insertions
        if isgap(refi)
            iszero(n_ins) && (markpos = refpos + 1)
            push!(insertion_buffer, seqi)
            n_ins += 1
        elseif !iszero(n_ins)
            seq = T(insertion_buffer)
            push!(edits, Edit{T,E}(Insertion(seq), UInt(markpos)))
            empty!(insertion_buffer)
            n_ins = 0
        end

        # Substitutions
        if !isgap(refi) && !isgap(seqi) && seqi != refi
            push!(edits, Edit{T,E}(Substitution{E}(seqi), UInt(refpos)))
        end
    end

    # Check for clips at the end of the alignment
    last_anchors = aln.a.aln.anchors[(end - 1):end]

    # Final indel, if applicable
    if !any(anchor -> anchor.op == OP_SOFT_CLIP, last_anchors)
        if !iszero(n_gaps)
            push!(edits, Edit{T,E}(Deletion(UInt(n_gaps)), UInt(markpos)))
        elseif !iszero(n_ins)
            push!(edits, Edit{T,E}(Insertion(T(insertion_buffer)), UInt(markpos)))
        end
    end

    return Haplotype(ref, edits)
end

"""
    _edits(h::Haplotype)

Gets the [`Edit`](@ref)s that comprise `h`
"""
_edits(h::Haplotype) = h.edits

"""
    reference(h::Haplotype)

Gets the reference sequence of `h`.
"""
reference(h::Haplotype) = h.ref
Base.:(==)(x::Haplotype, y::Haplotype) = x.ref == y.ref && x.edits == y.edits

"""
    reconstruct(h::Haplotype)

Apply the edits in `h` to the reference sequence of `h` and return the mutated sequence
"""
function reconstruct(h::Haplotype)
    len = length(reference(h)) + sum(edit -> _lendiff(edit), _edits(h))
    seq = copy(reference(h))
    resize!(seq, len % UInt)
    refpos = seqpos = 1
    for edit in _edits(h)
        while refpos < leftposition(edit)
            seq[seqpos] = reference(h)[refpos]
            refpos += 1
            seqpos += 1
        end
        editx = _mutation(edit)
        if editx isa Substitution
            seq[seqpos] = editx.x
            seqpos += 1
            refpos += 1
        elseif editx isa Deletion
            refpos += editx.len
        elseif editx isa Insertion
            for i in editx.seq
                seq[seqpos] = i
                seqpos += 1
            end
        end
    end
    while seqpos ≤ length(seq)
        seq[seqpos] = reference(h)[refpos]
        refpos += 1
        seqpos += 1
    end
    return seq
end
