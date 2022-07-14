# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Mutation types (`Substitution`, `Deletion`, `Insertion`)
- `Variant` type to store groups of mutations together
- `Variation` type to store and compare individual mutations
- `reconstruct!` function to build mutated sequences from `Variant`s
- `Variant` constructor to automatically detect mutations from a `BioAlignments.PairwiseAlignment`
- Methods to get reference and alternate bases from a `Variation`

[unreleased]: https://github.com/BioJulia/SequenceVariation.jl
