"""
Needs to be able to:
* Given a sequence and a reference, create a `Variant` that unambiguously represents
the sequence

* Given a `Variant` and a new reference, translate the variant to the new reference.

* Given a mutation and a reference and a sequence, determine if the sequence has that
mutation

TODO now:
* Create a string repr and parser for Edit, perhaps
    * A243T for sub
    * 119TAGGCTA for insertion
    * TGAGCTA9 for deletion
* Create a parser + print/show for edit
* Play around with some NGS results rel. to picked reference.
    * Is it easy to construct ref and variants? I.e. is API nice?
    * Is it nice and easy to check if a mut is present?
    *

* Implement "reference switching".
* Add tests
"""

using BioSequences
using BioAlignments
using SequenceVariation

const DNA_MODEL = BioAlignments.AffineGapScoreModel(EDNAFULL, gap_open=-25, gap_extend=-2)

align(a::BioSequence, b::BioSequence) = pairalign(GlobalAlignment(), a, b, DNA_MODEL).aln
seq1 = ungap!(dna"--ATGCGTGTTAGCAAC--TTATCGCG")
seq2 = ungap!(dna"TGATGCGTGT-AGCAACACTTATAGCG")
var = Variant(align(seq1, seq2))

@testset "VariationPosition" begin
    refseq = dna"ACAACTTTATCT"
    mutseq = dna"ACATCTTTATCT"

    read01 = AlignedSequence(mutseq[1:10], Alignment("10M", 1, 1))
    read02 = AlignedSequence(mutseq[3:12], Alignment("10M", 1, 3))
    aln01 = PairwiseAlignment(read01, refseq)
    aln02 = PairwiseAlignment(read02, refseq)
    @test Variant(aln01).edits == Variant(aln02).edits
end
