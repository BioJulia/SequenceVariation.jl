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
include("Variation.jl")

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
