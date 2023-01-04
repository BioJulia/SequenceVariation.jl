module SequenceVariation

"""
Needs to be able to:
* Given a sequence and a reference, create a `Haplotype` that unambiguously represents
the sequence

* Given a `Haplotype` and a new reference, translate the variant to the new reference.

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

using BioAlignments: BioAlignments, PairwiseAlignment, OP_SOFT_CLIP, sequence
using BioGenerics: BioGenerics, leftposition, rightposition
using BioSequences: BioSequences, BioSequence, NucleotideSeq, LongSequence, isgap
using BioSymbols: BioSymbol

export Deletion
export Haplotype
export Insertion
export Substitution
export Variation
export altbases
export mutation
export reconstruct!
export refbases
export reference
export translate
export variations

const BA = BioAlignments
const BS = BioSequences

struct Unsafe end
struct Inapplicable end

include("Edit.jl")
include("Haplotype.jl")
include("Variation.jl")

end # module
