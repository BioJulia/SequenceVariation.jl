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

## Variants

```@docs
Haplotype
reference(::Haplotype)
variations
reconstruct
```

## Variations

```@docs
Variation
reference(::Variation)
mutation
translate
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
