module SequenceVariation

"""
Needs to be able to:
* Given a sequence and a reference, create a `Variant` that unambiguously represents
the sequence

* Given a `Variant` and a new reference, translate the variant to the new reference.

"""

import BioSymbols: BioSymbol
using BioAlignments
using BioSequences

#=
import Automa
import Automa.RegExp: @re_str
=#

struct Unsafe end

const DEFAULT_AA_ALN_MODEL = AffineGapScoreModel(BLOSUM62, gap_open=-12, gap_extend=-2)
const DEFAULT_DNA_ALN_MODEL = AffineGapScoreModel(EDNAFULL, gap_open=-12, gap_extend=-2)
model(::AminoAcidSeq) = DEFAULT_AA_ALN_MODEL
model(::NucleotideSeq) = DEFAULT_DNA_ALN_MODEL

#=
const CONTEXT = Automa.CodeGenContext(
    vars=Automa.Variables(:p, :p_end, :p_eof, :ts, :te, :cs, :data, :mem, :byte),
    generator=:goto,
    checkbounds=false
)
=#

const BYTES = Union{String, SubString{String}, Vector{UInt8}}

"""
    Substitution

Represents the presence of a `T` at a given position. The position is stored
outside this struct.
"""
struct Substitution{T <: BioSymbol}
    x::T
end
Substitution(x::BioSymbol) = Substitution{typeof(x)}(x)

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

"""
    Insertion{S <: BioSequence}

Represents the insertion of a `S` into a sequence. The location of the insertion
is stored outside the struct.
"""
struct Insertion{S <: BioSequence}
    x::S

    function Insertion{S}(x::S) where {S <: BioSequence}
        isempty(x) && error("Insertion must be at least 1 symbol")
        new(x)
    end
end
Insertion(s::BioSequence) = Insertion{typeof(s)}(s)
Base.:(==)(x::Insertion, y::Insertion) = x.seq == y.seq

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

# Validate:
# A sequence is invalid if any of its operations are out of bounds, or the same position
# is affected by multiple edits.
function is_valid(v::Variant)
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
            pos in (first(valid_positions)+last_was_insert:last(valid_positions)) || return false
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

#=
# Substituion can only occur in 1:len
# Deletions of len L only in 1:(len-L)+1
# Insertions at 1:len+1
"If a Variant passes this validation, then it should be possible to unambiguously
reconstruct the variant."
function validate(v::Variant)
    len = length(v.ref)
    sort!(v.edits, by=x -> x.pos)
    lower = 1

    # No not allow two inserts right after each other at same pos, because then it's
    # ambiguous which to pick first.
    last_was_insert = false
    for edit in v.edits
        editx = edit.x
        pos = edit.pos
        if editx isa Substitution
            (pos < lower || pos > len) && error()
            lower = pos + 1
            last_was_insert = false
        elseif editx isa Deletion
            (pos < lower || pos > (len - editx.len + 1)) && error()
            lower = pos + editx.len
            last_was_insert = false
        else # Insertion
            (pos < lower + last_was_insert || pos > len + 1) && error()
            lower = pos
            last_was_insert = true
        end
    end
    return v
end
=#

# TODO: We NEED to include diffs at ends. Yes, they may not be proper mutations
# but we need it to reconstruct the actual seqs
function variant(seq::T, ref::T) where {T <: LongSequence{<:Union{AminoAcidAlphabet, NucleicAcidAlphabet}}}
    return variant(pairalign(OverlapAlignment(), seq, ref, model(seq)).aln)
end

function variant(aln::PairwiseAlignment{T, T}) where {T <: LongSequence{<:Union{AminoAcidAlphabet, NucleicAcidAlphabet}}}
    ref = aln.b
    E = eltype(typeof(ref))
    edits = Edit{T, E}[]
    result = Variant(ref, edits)
    refpos = seqpos = 0
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
        else
            if !iszero(n_gaps)
                push!(edits, Edit{T, E}(Deletion(UInt(n_gaps)), UInt(markpos)))
                n_gaps = 0
            end
        end

        # Check for insertions
        if isgap(refi)
            iszero(n_ins) && (markpos = refpos + 1)
            push!(insertion_buffer, seqi)
            n_ins += 1
        else
            if !iszero(n_ins)
                seq = T(insertion_buffer)
                push!(edits, Edit{T, E}(Insertion(seq), UInt(markpos)))
                empty!(insertion_buffer)
                n_ins = 0
            end
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
end

function Base.in(v::Variation, var::Variant)
    if v.ref != var.ref
        error("References must be equal")
    end
    any(v.edit == edit for edit in var.edits)
end

# How to check if a Variation is present in a sequence?
function hasvariation(seq::S, var::Variation{S, T}) where {S, T}
    var in variant(seq, var.ref)
end

end # module
