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
    valid, message = _is_valid(result)
    valid || error(message)
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
    _is_valid(h::Haplotype{S,T}) where {S,T}

Validate `h`. `h` is invalid if any of its operations are out of bounds, or the same
position is affected by multiple edits.
"""
function _is_valid(h::Haplotype{S,T}) where {S,T}
    # Empty references mean we have nothing to compare to
    isempty(reference(h)) && return (false, "Empty reference")
    # Empty edits simple means that this is the reference haplotype
    isempty(_edits(h)) && return (true, "")

    # It is valid for edits to exist within the space between the first edit and the entire
    # length of the reference. The reason we can't use simply use the length of the reference
    # is because we couldn't tell if a starting insertion overlapped with a future edit that
    # way. Also, force the range to be a signed integer to avoid subtraction overflow errors
    # when there is a starting insertion.
    valid_positions = UnitRange{Int128}(
        leftposition(first(_edits(h))), UInt64(length(reference(h)))
    )
    # Placeholder flag for iterating through multiple insertions at the same position
    last_was_insert = false

    for edit in _edits(h)
        pos = leftposition(edit)
        op = _mutation(edit)

        # Sanity check: for this to be a valid variant, it must be comprised of valid
        # variations
        _is_valid(Variation{S,T}(h.ref, edit, Unsafe())) ||
            return (false, "Invalid Variation")

        if op isa Substitution
            # For substitutions we simply do not allow another modification of the same base
            pos in valid_positions ||
                return (false, "Multiple modifications at same position")

            # Invalidate this base's and all previous positions. This only works because
            # we're working with a pre-sorted list of edits.
            valid_positions = (pos + 1):last(valid_positions)
            last_was_insert = false
        elseif op isa Insertion
            # Insertions affect 0 reference bases, so it does not modify the valid positions
            # for next op. However, we cannot have two insertions at the same position, because
            # then the order of them is ambiguous
            pos in (
                (first(valid_positions) - 1 + 2 * last_was_insert):(last(valid_positions) + 1)
            ) || return (false, "Multiple insertions at same position")

            last_was_insert = true
        elseif op isa Deletion
            # Deletions obviously invalidate the reference bases that are deleted.
            len = length(op)
            pos in (first(valid_positions):(last(valid_positions) - len + 1)) ||
                return (false, "Deletion out of range")

            valid_positions = (first(valid_positions) + len):last(valid_positions)
            last_was_insert = false
        end
    end
    return (true, "")
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

"""
    cigar(hap::Haplotype{S,T}) where {S,T}

Constructs a CIGAR string representing the alignment of the sequence of `hap` to its
reference.
"""
function BioAlignments.cigar(hap::Haplotype{S,T}) where {S,T}
    cigar_string = String[]

    mismatch_vars = filter(var -> !isa(mutation(var), Substitution), variations(hap))

    length(mismatch_vars) > 0 || return "$(length(reference(hap)))M"

    lastvar = first(mismatch_vars)

    leftposition(lastvar) > 1 && push!(cigar_string, "$(leftposition(lastvar))M")

    for var in mismatch_vars
        push!(cigar_string, _cigar_between(lastvar, var))
        push!(cigar_string, _cigar(var))
        lastvar = var
    end #for

    remaining_bases = length(reference(hap)) - rightposition(lastvar)
    remaining_bases > 0 && push!(cigar_string, "$(remaining_bases)M")

    return join(cigar_string, "")
end

"""
    translate(hap::Haplotype{S,T}, aln::PairwiseAlignment{S,S}) where {S,T}

Convert the variations in `hap` to a new reference sequence based upon `aln`. The alignment
rules follow the conventions of
[`translate(::Variation, PairwiseAlignment)`](@ref translate(::Variation{S,T}, ::PairwiseAlignment{S,S}) where {S,T}).
Indels at the beginning or end may not be preserved. Returns a new
[`Haplotype`](@ref)
"""
function translate(hap::Haplotype{S,T}, aln::PairwiseAlignment{S,S}) where {S,T}
    vars = variations(hap)
    new_ref = BA.sequence(aln)
    translated_vars = Variation{S,T}[]
    for v in vars
        new_v = translate(v, aln)
        isnothing(new_v) || push!(translated_vars, new_v)
    end
    return Haplotype(new_ref, translated_vars)
end
