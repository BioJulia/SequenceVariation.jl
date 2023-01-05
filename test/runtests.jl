"""
Needs to be able to:
* Given a sequence and a reference, create a `Haplotype` that unambiguously represents
the sequence

* Given a `Haplotype` and a new reference, translate the variant to the new reference.

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

using Aqua
using BioAlignments:
    AffineGapScoreModel,
    GlobalAlignment,
    EDNAFULL,
    pairalign,
    Alignment,
    AlignedSequence,
    PairwiseAlignment
using BioSequences: BioSequence, @dna_str, ungap!
using BioSymbols: DNA_A
using SequenceVariation
using Test

const DNA_MODEL = AffineGapScoreModel(EDNAFULL; gap_open=-25, gap_extend=-2)

align(a::BioSequence, b::BioSequence) = pairalign(GlobalAlignment(), a, b, DNA_MODEL).aln
seq1 = ungap!(dna"--ATGCGTGTTAGCAAC--TTATCGCG")
seq2 = ungap!(dna"TGATGCGTGT-AGCAACACTTATAGCG")
seq3 = ungap!(dna"-GATGCGTGT-AGCAACACTTATCGC-")
var = Haplotype(align(seq1, seq2))

@testset "EditSorting" begin
    S = typeof(seq1)
    T = eltype(seq1)
    @test SequenceVariation.Edit{S,T}(Deletion(1), 1) <
        SequenceVariation.Edit{S,T}(Deletion(1), 2)
    @test SequenceVariation.Edit{S,T}(Deletion(1), 1) <
        SequenceVariation.Edit{S,T}(Deletion(2), 1)
end

@testset "HaplotypeRoundtrip" begin
    for v in variations(var)
        @test v in var
        @test v in Haplotype(seq2, [v])
    end
end

@testset "HaplotypeReconstruction" begin
    @test reconstruct(var) == seq1
end

@testset "VariationSorting" begin
    @test Variation(seq2, "A3T") < Variation(seq2, "T4A")
end

@testset "HaplotypeTranslation" begin
    ref1 = seq2
    ref2 = seq3
    alt = seq1

    @test all(
        v -> v in Haplotype(align(alt, ref2)), variations(translate(var, align(ref2, ref1)))
    )
end

@testset "VariationPosition" begin
    refseq = dna"ACAACTTTATCT"
    mutseq = dna"ACATCTTTATCT"

    read01 = AlignedSequence(mutseq[1:10], Alignment("10M", 1, 1))
    read02 = AlignedSequence(mutseq[3:12], Alignment("10M", 1, 3))
    aln01 = PairwiseAlignment(read01, refseq)
    aln02 = PairwiseAlignment(read02, refseq)
    @test Haplotype(aln01).edits == Haplotype(aln02).edits
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
    var = Haplotype(aln)

    sub = Variation(refseq, "A4T")
    @test first(variations(var)) == sub
end

@testset "VariationTranslation" begin
    ref1 = seq2
    ref2 = seq3
    ref3 = copy(ref1)
    ref3[1] = DNA_A
    alt = seq1

    ref2_on_ref1 = align(ref2, ref1)
    alt_on_ref1 = align(alt, ref1)
    ref1_on_alt = align(ref1, alt)
    ref3_on_ref1 = align(ref3, ref1)

    # Test insertion at position zero
    @test translate(Variation(ref1, "0CAT"), ref2_on_ref1) ==
        SequenceVariation.Inapplicable()
    @test translate(Variation(ref1, "0CAT"), ref3_on_ref1) == Variation(ref3, "0CAT")

    # Test substitution on deletion
    @test isnothing(translate(Variation(ref1, "A17T"), ref1_on_alt))

    # Test substitution to itself
    @test isnothing(translate(Variation(ref1, "A23C"), ref2_on_ref1))

    # Test simple substitution
    @test translate(Variation(ref1, "T10A"), ref2_on_ref1) == Variation(ref2, "T9A")

    # Test simple deletion
    @test translate(Variation(ref1, "Δ17-18"), ref2_on_ref1) == Variation(ref2, "Δ16-17")

    # Test deletion on deletion
    @test isnothing(translate(Variation(ref1, "Δ17-17"), alt_on_ref1))

    # Test simple insertion
    @test translate(Variation(ref1, "6CAT"), ref2_on_ref1) == Variation(ref2, "5CAT")

    # Test two insertions at the same position
    @test translate(Variation(ref1, "10G"), alt_on_ref1) == SequenceVariation.Inapplicable()

    # Test insertion on deletion
    @test translate(Variation(ref1, "17CAT"), alt_on_ref1) ==
        SequenceVariation.Inapplicable()
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

@testset "SoftclipHaplotype" begin
    refseq = dna"GATTACA"
    mutseq = dna"GATTACAAAA"

    refvar = Haplotype(refseq, SequenceVariation.Edit{typeof(refseq),eltype(refseq)}[])

    # Test for ending soft clip
    @test Haplotype(
        PairwiseAlignment(AlignedSequence(mutseq, Alignment("7=3S", 1, 1)), refseq)
    ) == refvar

    # Test for ending soft+hard clip
    @test Haplotype(
        PairwiseAlignment(AlignedSequence(mutseq, Alignment("7=3S2H", 1, 1)), refseq)
    ) == refvar

    # Test that ending insertions are still valid
    @test length(
        Haplotype(
            PairwiseAlignment(AlignedSequence(mutseq, Alignment("7=3I", 1, 1)), refseq)
        ).edits,
    ) == 1

    # Test that out-of-bounds bases are still caught
    @test_throws BoundsError Haplotype(
        PairwiseAlignment(AlignedSequence(mutseq, Alignment("7=3X", 1, 1)), refseq)
    )
end

@testset "Aqua" begin
    Aqua.test_ambiguities(SequenceVariation; recursive=false)
    # TODO: Refactor `Edit` so that this test doesn't fail
    # TODO: This test _should_ be set to @test_fails, but Aqua's syntax doesn't allow that
    # Aqua.test_unbound_args(SequenceVariation)
    Aqua.test_undefined_exports(SequenceVariation)
    Aqua.test_piracy(SequenceVariation)
    Aqua.test_project_extras(SequenceVariation)
    Aqua.test_stale_deps(SequenceVariation)
    Aqua.test_deps_compat(SequenceVariation)
    Aqua.test_project_toml_formatting(SequenceVariation)
end
