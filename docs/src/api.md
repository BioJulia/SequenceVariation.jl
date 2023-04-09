```@meta
CurrentModule = SequenceVariation
DocTestSetup = quote
    using SequenceVariation
end
```

# API Reference

## Edits

```@docs
Substitution
Deletion
Insertion
```

## Haplotypes

```@docs
Haplotype
reference(::Haplotype)
variations
reconstruct
BioAlignments.alignment
BioAlignment.cigar
translate(::Haplotype{S,T}, ::PairwiseAlignment{S,S}) where {S,T}
```

## Variations

```@docs
Variation
reference(::Variation)
mutation
translate(::Variation{S,T}, ::PairwiseAlignment{S,S}) where {S,T}
refbases
altbases
```

## Private API

### Edits

```@docs
Edit
_mutation
_lendiff
```

### Variants

```@docs
_edits
_is_valid(::Haplotype)
```

### Variations

```@docs
_edit
_is_valid(::Variation)
```
