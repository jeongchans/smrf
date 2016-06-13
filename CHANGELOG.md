# Change Log
All notable changes to this project will be documented in this file.


## Unreleased

### Changed
- Use PMRF v0.4.0, much faster and versatile.
- MRF file is written in PMRF v0.4.0 compatible format.

### Compatibility Warning
- Not support PMRF v0.2.x.


## 0.2.2 - 2016-06-02

### Added
- Specify the edge file for MRF architecture. Use `--edge` option.
- Specify the mrf file for SMRF model. Use `--mrf` option.
- Show the release version with `--version` option

### Changed
- Change the command-line options. Input files are required, but output files are optional.
- Require <msa_file> and <pdb_file>.
- Specify the output file for the positional coevolution scores. Use `--pos` option. If not specified, stdout will be used as default.
- Specify the output file for the pairwise coevolution scores. Use `--pair` option. If not specified, the pairwise coevoltion scores will not be reported.
- Show help message with `-h` and `--help` options.


## 0.2.1 - 2016-02-04

### Added
- This CHANGELOG file to contain notable changes.

### Changed
- Organize the installation procedure in README document.

### Fixed
- Fix the serious bug in VERSION assignment.
- Fix the initial generation of configuration file.
- Fix the link for the PMRF software in the README document.


## 0.2.0 - 2016-01-07 [YANKED]

### Added
- Build the structure-based Markov random field.
- Calculate the pairwise coevolution scores.
- Calculate the positional coevolution scores.
^