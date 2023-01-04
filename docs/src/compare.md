```@meta
CurrentModule = SequenceVariation
```

# Comparing variations in sequences

## Checking for variations in a known variant

Looking for a known [`Variation`](@ref) within a [`Haplotype`](@ref) is
efficiently accomplished using the `in` operator.

```@setup call_variants
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

```@example call_variants
println("\tOvis aires\tHomo sapiens")
for v in vcat(variations(bos_ovis_variant), variations(bos_human_variant))
    is_sheep = v in bos_ovis_variant
    is_human = v in bos_human_variant
    println("$v\t$is_sheep\t\t$is_human")
end
```

## Constructing new variants based on other variations

New variants can be constructed using variations. This might be useful to pool
variations found on different reads or to filter variations from a variant
that aren't validated by another variant.

```@repl call_variants
sheeple = vcat(variations(bos_ovis_variant), variations(bos_human_variant));
Haplotype(bovine, sheeple)
reconstruct!(bovine, ans)
```
