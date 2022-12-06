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

using BioAlignments
using BioSequences
using SequenceVariation
using Test

const DNA_MODEL = BioAlignments.AffineGapScoreModel(EDNAFULL, gap_open=-25, gap_extend=-2)

align(a::BioSequence, b::BioSequence) = pairalign(GlobalAlignment(), a, b, DNA_MODEL).aln
seq1 = ungap!(dna"--ATGCGTGTTAGCAAC--TTATCGCG")
seq2 = ungap!(dna"TGATGCGTGT-AGCAACACTTATAGCG")
var = Variant(align(seq1, seq2))

@testset "VariantRoundtrip" begin
    for v in variations(var)
        @test v in var
        @test v in Variant(seq2, [v])
    end
end

@testset "VariationPosition" begin
    refseq = dna"ACAACTTTATCT"
    mutseq = dna"ACATCTTTATCT"

    read01 = AlignedSequence(mutseq[1:10], Alignment("10M", 1, 1))
    read02 = AlignedSequence(mutseq[3:12], Alignment("10M", 1, 3))
    aln01 = PairwiseAlignment(read01, refseq)
    aln02 = PairwiseAlignment(read02, refseq)
    @test Variant(aln01).edits == Variant(aln02).edits
end

@testset "VariationParsing" begin
    refseq = dna"ACAACTTTATCT"

    sub = Variation(refseq, "A4T")
    del = Variation(refseq, "Δ4-5")
    ins = Variation(refseq, "4TT")

    @test mutation(sub) isa Substitution
    @test mutation(del) isa Deletion
    @test mutation(ins) isa Insertion
end

@testset "VariationRetrieval" begin
    refseq = dna"ACAACTTTATCT"
    mutseq = dna"ACATCTTTATCT"

    read = AlignedSequence(mutseq[1:10], Alignment("10M", 1, 1))
    aln = PairwiseAlignment(read, refseq)
    var = Variant(aln)

    sub = Variation(refseq, "A4T")
    @test first(variations(var)) == sub
end

@testset "VariationBases" begin
    # Test substition bases
    @test refbases(Variation(dna"ATCGA", "C3G")) == dna"C"
    @test altbases(Variation(dna"ATCGA", "C3G")) == dna"G"

    # Test single deletion bases
    @test refbases(Variation(dna"ATCGA", "Δ3-3")) == dna"TC"
    @test altbases(Variation(dna"ATCGA", "Δ3-3")) == dna"T"

    # Test multiple deletion bases
    @test refbases(Variation(dna"ATCGA", "Δ3-4")) == dna"TCG"
    @test altbases(Variation(dna"ATCGA", "Δ3-4")) == dna"T"

    # Test first position deletion
    @test refbases(Variation(dna"ATCGA", "Δ1-1")) == dna"AT"
    @test altbases(Variation(dna"ATCGA", "Δ1-1")) == dna"T"

    # Test single insertion bases
    @test refbases(Variation(dna"ATCGA", "3A")) == dna"C"
    @test altbases(Variation(dna"ATCGA", "3A")) == dna"CA"

    # Test multiple insertion bases
    @test refbases(Variation(dna"ATCGA", "3TAG")) == dna"C"
    @test altbases(Variation(dna"ATCGA", "3TAG")) == dna"CTAG"

    # Test first position insertion
    @test refbases(Variation(dna"ATCGA", "1C")) == dna"A"
    @test altbases(Variation(dna"ATCGA", "1C")) == dna"CA"
end

@testset "SoftclipVariant" begin
    refseq = dna"GATTACA"
    mutseq = dna"GATTACAAAA"

    refvar = Variant(refseq, SequenceVariation.Edit{typeof(refseq), eltype(refseq)}[])

    # Test for ending soft clip
    @test Variant(PairwiseAlignment(AlignedSequence(mutseq, Alignment("7=3S", 1, 1)), refseq)) == refvar

    # Test for ending soft+hard clip
    @test Variant(PairwiseAlignment(AlignedSequence(mutseq, Alignment("7=3S2H", 1, 1)), refseq)) == refvar
end
