```@meta
CurrentModule = SequenceVariation
```

# Working with individual variations

## Construction

Individual [`Variation`](@ref)s can be made using a reference sequence and
string syntax

| Variation type | Syntax | Interpretation | Example |
|:--- |:--- |:--- |:--- |
| Substitutions | `<REF><POS><ALT>` | `<ALT>` is substituted for `<REF>` in position `<POS>` | `"G16C"` |
| Deletions | `Δ<START>-<END>` | All bases (inclusive) between `<START>` and `<END>` are deleted. It is valid to have `<START>` equal `<END>`: that is a deletion of one base. | `"Δ1-2"` |
| Insertions | `<POS><ALT>` | `<ALT>` is inserted between positions `<POS>` and `<POS>+1` | `"11T"` |

```@repl
using BioSequences: @dna_str
using SequenceVariation
bovine_ins = dna"GACCGGCTGCATTCGAGGCTGCCAGCAAGCAG"
Variation(bovine_ins, "C4A")
mutation(ans)
typeof(mutation(Variation(bovine_ins, "Δ13-14")))
mutation(Variation(bovine_ins, "25ACA"))
```

## Extraction

Sequence variations may also be extracted wholesale from a [`Variant`](@ref)
using the [`variations`](@ref) function.

```@setup call_variants
using SequenceVariation, BioAlignments, BioSequences

bovine = dna"GACCGGCTGCATTCGAGGCTGCCAGCAAGCAG";
ovine  = dna"GACCGGCTGCATTCGAGGCTGTCAGCAAACAG";
human  = dna"GACAGGCTGCATCAGAAGAGGCCATCAAGCAG";

bos_ovis_alignment =
    PairwiseAlignment(AlignedSequence(ovine, Alignment("32M", 1, 1)), bovine);
bos_human_alignment =
    PairwiseAlignment(AlignedSequence(human, Alignment("32M", 1, 1)), bovine);

bos_ovis_variant = Variant(bos_ovis_alignment)
bos_human_variant = Variant(bos_human_alignment)
```

```@repl call_variants
variations(bos_ovis_variant)
variations(bos_human_variant)
```

## Reference switching

An individual variation can be mapped to a new reference sequence given an
alignment between the new and old references using the [`translate`](@ref).

```@repl call_variants
ovis_human_alignment =
    PairwiseAlignment(AlignedSequence(human, Alignment("32M", 1, 1)), ovine)
human_variation = first(variations(bos_ovis_variant))
reference(ans) == bovine
SequenceVariation.translate(human_variation, ovis_human_alignment)
reference(ans) == bovine
```
