```@meta
CurrentModule = SequenceVariation
```

# Comparing variations in sequences

## Checking for variations in a known haplotype

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

bos_ovis_haplotype = Haplotype(bos_ovis_alignment)
bos_human_haplotype = Haplotype(bos_human_alignment)
```

```@example call_variants
println("\tOvis aires\tHomo sapiens")
for v in vcat(variations(bos_ovis_haplotype), variations(bos_human_haplotype))
    is_sheep = v in bos_ovis_haplotype
    is_human = v in bos_human_haplotype
    println("$v\t$is_sheep\t\t$is_human")
end
```

## Constructing new haplotypes based on other variations

New haplotypes can be constructed using variations. This might be useful to pool
variations found on different reads or to filter variations from a haplotype
that aren't validated by another haplotype.

```@repl call_variants
sheeple = vcat(variations(bos_ovis_haplotype), variations(bos_human_haplotype));
Haplotype(bovine, sheeple)
reconstruct!(bovine, ans)
```
