
# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2025-08-11

### Added
- R analysis container (r-analysis:1.0) with R 4.3.2 and required packages for pipeline visualization
- Support for BCF output format from Minimac4 imputation
- Test configuration for v6 phased data (v6_chr21_phased.config)

### Fixed
- Minimac4 output format issue - changed from SAV to BCF format for better compatibility
- Reference panel chromosome placeholder replacement in generate_frequency process
- Chromosome naming issues with double "chr" prefix in reference panels
- Minimac4 parameters for MSAV format compatibility

### Changed
- Updated impute_minimac4 process to output BCF format instead of SAV
- Modified generate_frequency workflow to properly handle chromosome-specific reference panel paths
- Adjusted minRatio parameter to 0.01 for better variant ratio handling

## [0.1] - 2016-08-22
Initial release of h3achipimputation pipeline, created as part of the  H3ABioNet Hackathon held in Pretoria, SA in 2016.
