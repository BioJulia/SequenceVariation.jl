# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.2] - 2023-01-28

### Fixed

- `Haplotype` validation ([#36](https://github.com/BioJulia/SequenceVariation.jl/issues/36)/[#37](https://github.com/BioJulia/SequenceVariation.jl/pull/37))

## [0.2.1] - 2023-01-11

### Added

- `translate` functionality for `Haplotype`s ([#31](https://github.com/BioJulia/SequenceVariation.jl/pull/31))
- More informative errors when constructing invalid `Haplotype`s ([#34](https://github.com/BioJulia/SequenceVariation.jl/pull/34))

## [0.2.0] - 2023-01-10

### Added

- Tutorial-type documentation ([#28](https://github.com/BioJulia/SequenceVariation.jl/pull/28))
- `Base.isless` implementation for `Edit` and `Variation` ([#32](https://github.com/BioJulia/SequenceVariation.jl/pull/32))

### Changed

- Code now follows [Blue style](https://github.com/invenia/BlueStyle) ([#28](https://github.com/BioJulia/SequenceVariation.jl/pull/28))
- :bomb: [BREAKING] Public and private API defined based on Blue style guidelines ([#28](https://github.com/BioJulia/SequenceVariation.jl/pull/28))
- :bomb: [BREAKING] Renamed type `Variant` to `Haplotype` ([#20](https://github.com/BioJulia/SequenceVariation.jl/issues/20)/[#29](https://github.com/BioJulia/SequenceVariation.jl/pull/29))
- :bomb: [BREAKING] Refactored `reconstruct` function to use `Haplotype`'s reference sequence ([#30](https://github.com/BioJulia/SequenceVariation.jl/pull/30))

### Removed

- Windows and MacOS CI tests ([#28](https://github.com/BioJulia/SequenceVariation.jl/pull/28))

## [0.1.4] - 2022-12-17

### Fixed

- Soft clips at end of alignment causing invalid `Variant`s ([#25](https://github.com/BioJulia/SequenceVariation.jl/issues/25)/[#26](https://github.com/BioJulia/SequenceVariation.jl/pull/26))

## [0.1.3] - 2022-11-22

### Changed

- Variations getter now returns type-parameterized vector ([#23](https://github.com/BioJulia/SequenceVariation.jl/pull/23))

## [0.1.2] - 2022-10-04

### Changed

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

[unreleased]: https://github.com/BioJulia/SequenceVariation.jl/compare/v0.2.2...HEAD
[0.2.2]: https://github.com/BioJulia/SequenceVariation.jl/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/BioJulia/SequenceVariation.jl/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/BioJulia/SequenceVariation.jl/compare/v0.1.4...v0.2.0
[0.1.4]: https://github.com/BioJulia/SequenceVariation.jl/compare/v0.1.3...v0.1.4
[0.1.3]: https://github.com/BioJulia/SequenceVariation.jl/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/BioJulia/SequenceVariation.jl/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/BioJulia/SequenceVariation.jl/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/BioJulia/SequenceVariation.jl/releases/tag/v0.1.0
