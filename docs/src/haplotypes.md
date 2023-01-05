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
human2 = reconstruct(bos_human_haplotype)
human2 == bovine
human2 == human
```

## Reference switching

All variations within a haplotype can be mapped to a new reference sequence
given an alignment between the new and old references using the
[`translate`](@ref translate(::Haplotype{S,T}, ::PairwiseAlignment{S,S}) where {S,T})
function. This could be useful if variants were called against a reference
sequence for the entire species, but need to be analyzed as variants of a
subtype later.

```@repl call_variants
ovis_human_alignment =
    PairwiseAlignment(AlignedSequence(human, Alignment("32M", 1, 1)), ovine)
SequenceVariation.translate(bos_ovis_haplotype, ovis_human_alignment)
```
