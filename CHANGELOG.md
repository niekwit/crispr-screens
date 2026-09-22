# Changelog

## [1.0.0](https://www.github.com/niekwit/crispr-screens/compare/v0.10.0...v1.0.0) (2026-09-22)


### ⚠ BREAKING CHANGES

* Align trimmed reads to the sgRNA index with Bowtie (v1) instead of HISAT2. The count rule now uses the shell directive instead of count.sh, and the alignment rate plot parses the Bowtie log.

### Features

* add script to annotate sgRNAs with genomic coordinates for CRISPRcleanR ([6b0f905](https://www.github.com/niekwit/crispr-screens/commit/6b0f905a8d5ab635e706ec8cc5126c5218b29b69))
* redesign sgRNA rank plot with density overview and enriched/depleted split ([8d327ce](https://www.github.com/niekwit/crispr-screens/commit/8d327ceeec0619e1e22b7f3948f2b9ad31ecdf55))
* replace HISAT2 with Bowtie for sgRNA alignment ([6238da9](https://www.github.com/niekwit/crispr-screens/commit/6238da9ca0272740f07487662c3e3a0c06706120))


### Bug Fixes

* grant contents/pull-requests write permission to release-please ([b218d02](https://www.github.com/niekwit/crispr-screens/commit/b218d022ce12148d8795130d6d3d8e6196dc3ca0))
* handle genome-scope runs without gene column and clarify errors ([6c9106e](https://www.github.com/niekwit/crispr-screens/commit/6c9106e52dbbbc88381497c6bc6c128027869b0e))
* prevent bowtie crash on short reads and ggplot2 facet crash in count distribution plot ([189c162](https://www.github.com/niekwit/crispr-screens/commit/189c16238227e37597d6bd77956a67b10e4a89bd))
