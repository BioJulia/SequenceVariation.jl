# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

- ## Changed

- Updated dependency compats ([#21](https://github.com/BioJulia/SequenceVariation.jl/pull/21))
  - BioAlignments: 2 -> 2,3
  - BioSequences: 2 -> 2,3
  - BioSymbols: 4 -> 4,5

## [0.1.1] - 2022-07-21

### Added

- Constructor for `Variant` based on `Variation`s ([#18](https://github.com/BioJulia/SequenceVariation.jl/pull/18))

## [0.1.0] - 2022-07-14

### Added

- Mutation types (`Substitution`, `Deletion`, `Insertion`)
- `Variant` type to store groups of mutations together
- `Variation` type to store and compare individual mutations
- `reconstruct!` function to build mutated sequences from `Variant`s
- `Variant` constructor to automatically detect mutations from a `BioAlignments.PairwiseAlignment`
- Methods to get reference and alternate bases from a `Variation`

[unreleased]: https://github.com/BioJulia/SequenceVariation.jl/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/BioJulia/SequenceVariation.jl/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/BioJulia/SequenceVariation.jl/releases/tag/v0.1.0
