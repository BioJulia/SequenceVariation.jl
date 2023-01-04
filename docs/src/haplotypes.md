```@meta
CurrentModule = SequenceVariation
```

# Working with variants

## Calling variants

The first step in working with sequence variation is to identify (call)
variations. SequenceVariation can directly call variants using the
`Haplotype(::PairwiseAlignment)` constructor of the [`Haplotype`](@ref) type.

```@repl call_variants
using SequenceVariation, BioAlignments, BioSequences

bovine = dna"GACCGGCTGCATTCGAGGCTGCCAGCAAGCAG";
ovine  = dna"GACCGGCTGCATTCGAGGCTGTCAGCAAACAG";
human  = dna"GACAGGCTGCATCAGAAGAGGCCATCAAGCAG";

bos_ovis_alignment =
    PairwiseAlignment(AlignedSequence(ovine, Alignment("32M", 1, 1)), bovine);
bos_human_alignment =
    PairwiseAlignment(AlignedSequence(human, Alignment("32M", 1, 1)), bovine);

bos_ovis_variant = Haplotype(bos_ovis_alignment)
bos_human_variant = Haplotype(bos_human_alignment)
```

## Sequence reconstruction

If the alternate sequence of a variant is no longer available (as is often the
case when calling variants from alignment files), then the sequence can be
retrieved using the [`reconstruct!`](@ref) function.

```@repl call_variants
human2 = copy(bovine);
reconstruct!(human2, bos_human_variant)
human2 == bovine
human2 == human
```
