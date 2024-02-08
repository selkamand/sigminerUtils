
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sigminerUtils

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/sigminerUtils)](https://CRAN.R-project.org/package=sigminerUtils)
<!-- badges: end -->

Utility package to simplify the process of running sigminer, storing
signature results and producing key visualisations

*Warning:* This package is in early development and not ready for use

## Installation

You can install the development version of sigminerUtils like so:

``` r
# install.packages('remotes')
remotes::install_github('selkamand/sigminerUtils')
```

## Usage

1.  Start by creating a Signature database

``` r
sig_create_database('./my_sig_database.sqlite', overwrite = TRUE)
```

2.  Run Signature Analysis

``` r
path_maf <- system.file(package = "sigminerUtils", path = "test.maf")

sig_analyse_mutations(maf = path_maf, ref = "hg19")
```

\`\`

3.  Load results from output folder into the database

``` r
sig_add_to_database('./signatures', './my_sig_database.sqlite', ref = "hg19")
```

Note you can optionally supply sample metadata

``` r
sample_metadata <- data.frame(
  sampleId = c('TCGA-A2-A0T5-01', 'TCGA-CA-6717-01', 'TCGA-CF-A9FF-01'),
  disease = c('Breast', 'Colorectal', 'Bladder'),
  description = character(3)
)

sig_add_to_database('./signatures', './my_sig_database.sqlite', ref = "hg19", metadata = sample_metadata)
```
