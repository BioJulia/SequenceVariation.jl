```@meta
CurrentModule = SequenceVariation
```

# Working with haplotypes

## Calling variants

The first step in working with sequence variation is to identify (call)
variations between two sequences. SequenceVariation can directly call variants
using the `Haplotype(::PairwiseAlignment)` constructor of the
[`Haplotype`](@ref) type.

```@repl call_variants
using SequenceVariation, BioAlignments, BioSequences

bovine = dna"GACCGGCTGCATTCGAGGCTGCCAGCAAGCAG";
ovine  = dna"GACCGGCTGCATTCGAGGCTGTCAGCAAACAG";
human  = dna"GACAGGCTGCATCAGAAGAGGCCATCAAGCAG";

bos_ovis_alignment =
    PairwiseAlignment(AlignedSequence(ovine, Alignment("32M", 1, 1)), bovine);
bos_human_alignment =
    PairwiseAlignment(AlignedSequence(human, Alignment("32M", 1, 1)), bovine);

bos_ovis_haplotype = Haplotype(bos_ovis_alignment)
bos_human_haplotype = Haplotype(bos_human_alignment)
```

## Sequence reconstruction

If the alternate sequence of a haplotype is no longer available (as is often the
case when calling variants from alignment files), then the sequence can be
retrieved using the [`reconstruct`](@ref) function.

```@repl call_variants
human2 = copy(bovine);
reconstruct(human2, bos_human_haplotype)
human2 == bovine
human2 == human
```
