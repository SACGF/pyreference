## Unreleased

### Changed

## [0.7.3] - 2023-07-06

### Changed

- Fix gene.representative_transcript dying with "AttributeError: module 'sys' has no attribute 'maxint'" in Python3

## [0.7.2] - 2022-11-21

### Added

- New Gene properties 'description', 'summary', 'map_location' and 'biotype'
- Support for Fasta reference genomes that use contigs for sequence names (eg NCBI)

### Changed

- We now use [cdot](https://github.com/SACGF/cdot) JSON.gz files

## [0.6.3] - 2022-01-12

### Changed

- Fixed bug where pyreference_biotype.py crashed due to args not matching method signature

## [0.6.2] - 2022-01-12

### Added

- Include pyreference_biotype.py script in PyPi distribution
- Removed individual graphs, improved appearance of stacked bar graph

### Changed

- Fixes for pyreference_biotype, pin HTSeq version to stop crash

## [0.6] - 2021-11-05

### Added

- Handle Ensembl specific GTFs
- Support for GFF3
- Store gene/transcript versions
- Store HGNC, description, cDNA_match (refseq transcript/genome alignment gaps)
- Store URL (where GTF/GFF was downloaded from eg RefSeq/Ensembl FTP site)

### Changed

- Fix for deprecated BioPython code

## [0.5] - 2020-02-24

### Changed

- Fixed Python 3.7 issue - ConfigParser mandatory arguments

## [0.4] - 2019-10-31

### Changed

- Fix Python3 issues
- Use PySam instead of PyFasta (performance issues at high chromosome coordinates)

## [0.3] - 2018-01-26

## Added

- Store GTF/GFF path and md5sum in JSON

### Changed

- TSS uses representative transcript start rather than most 3' transcript start

## [0.2] - 2018-01-25

### Added

- Be able to retrieve multiple genes at a time via list
- Option to decompress Gzip in memory to get around server shared filesystem issues

### Removed

- Removed non-standard chromosomes

## [0.1] - 2018-01-24

### Added

- Initial commit. Created project, extracted existing code from SACGF bioinformatics repo
- Wrote GTF to JSON converter and loader

[unreleased]: https://github.com/SACGF/pyreference/compare/v0.6.3...HEAD
[0.6.3]: https://github.com/SACGF/pyreference/compare/v0.6.2...v0.6.3
[0.6.2]: https://github.com/SACGF/pyreference/compare/v0.6...v0.6.2
[0.6]: https://github.com/SACGF/pyreference/compare/v0.5...v0.6
[0.5]: https://github.com/SACGF/pyreference/compare/v0.4...v0.5
[0.4]: https://github.com/SACGF/pyreference/compare/v0.3...v0.4
[0.3]: https://github.com/SACGF/pyreference/compare/v0.2...v0.3
[0.2]: https://github.com/SACGF/pyreference/compare/v0.1...v0.2
[0.1]: https://github.com/SACGF/pyreference/releases/tag/v0.1
