
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sigminerUtils

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/sigminerUtils)](https://CRAN.R-project.org/package=sigminerUtils)
[![Codecov test
coverage](https://codecov.io/gh/selkamand/sigminerUtils/branch/master/graph/badge.svg)](https://app.codecov.io/gh/selkamand/sigminerUtils?branch=master)
[![R-CMD-check](https://github.com/selkamand/sigminerUtils/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/selkamand/sigminerUtils/actions/workflows/R-CMD-check.yaml)
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

Load required packages

``` r
library(sigminerUtils)
```

1.  Start by creating a Signature database

``` r
sig_create_database('./my_sig_database.sqlite', overwrite = TRUE)
```

2.  Run Signature Analysis

``` r
# Variant Files
path_maf <-  system.file(package = "sigminerUtils", path = "test.maf") 
path_copynumber <- system.file(package = "sigminerUtils", path = "test.copynumber.csv") 
path_sv <- system.file(package = "sigminerUtils", path = "test.sv.csv") 

sig_analyse_mutations(
  maf = path_maf, 
  copynumber = read.csv(path_copynumber),
  structuralvariant = read.csv(path_sv),
  ref = "hg19", 
  output_dir = "./signatures"
  )
```

3.  Create a reference matrix

``` r
sig_create_reference_set(
  path_to_signature_directory = "./signatures",
  outfolder = "./reference"
  )
```

4.  Use reference matrix to rerun your signature analysis

Using a reference matrix allows samples to be compared to a reference
database

``` r
# Variant Files
path_maf <-  system.file(package = "sigminerUtils", path = "test.maf")
path_copynumber <- system.file(package = "sigminerUtils", path = "test.copynumber.csv") 
path_sv <- system.file(package = "sigminerUtils", path = "test.sv.csv") 

# Database of reference of sample catalogues
path_reference_tallies <- system.file(package = "sigminerUtils", path = "reference/refmatrix.tally.parquet")

sig_analyse_mutations(
  maf = path_maf, 
  copynumber = read.csv(path_copynumber),
  structuralvariant = read.csv(path_sv),
  ref = "hg19", 
  ref_tallies = path_reference_tallies,
  output_dir = "./signatures"
)
```

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

### Run in Single Sample Mode from Oncoanalyser outputs

``` r
path_snvs <- system.file("colo829_testfiles/COLO829v003T.purple.somatic.vcf.gz", package = "sigminerUtils")
path_cnvs <- system.file("colo829_testfiles/COLO829v003T.purple.cnv.somatic.tsv", package = "sigminerUtils")
path_svs <-  system.file("colo829_testfiles/COLO829v003T.purple.sv.vcf.gz", package = "sigminerUtils")

sig_analyse_mutations_single_sample_from_files(
  sample_id = "COLO829v003T", vcf_snv = path_snvs, segment = path_cnvs, vcf_sv = path_svs,
  pass_only = TRUE,
  ref = "hg38",
  output_dir = "colo829_signature_results"
)
```

### Run in Cohort Mode using a manifest file

Need to find a way to run this when manifest includes

``` r
# Extract PCAWG data into working directory
pcawg_tar = system.file("pcawg/pcawg_data.tar.gz", package = "sigminerUtils")
untar(tarfile = pcawg_tar, exdir = "pcawg")

# Supply a manifest file
path_manifest <- system.file("pcawg/manifest.tsv", package = "sigminerUtils")

# Run the signature analysis
sig_analyse_cohort_from_files(
  manifest = path_manifest,
  ref = "hg19",
  output_dir = "pcawg_signature_results", 
  cores=1
)
#> ✔ Found required columns in manifest: [sample, copynumber, and snv]
#> ℹ Running signature analysis for [5] samples. This may take a while ...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Found 1 samples described in the VCF [DO1000]
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ Returning data for only sample [DO1000] samples in format column of VCF.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Running Signature Analysis. This will take some time ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Mutational Signature Analysis ───────────────────────────────────────────────
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Checking arguments ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ Output Directory: 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ Reference Genome: BSgenome.Hsapiens.UCSC.hg19
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'Registered S3 method overwritten by 'sigminer':
#>   method      from
#>   print.bytes Rcpp
#> Please note that the generated MAF object is designed for mutational signature analysis, not recommended for Maftools analysis!
#> 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Tally (Small Variants) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:23.733146]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:26.874015]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> Attaching package: 'S4Vectors'
#> 
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> 
#> The following objects are masked from 'package:base':
#> 
#>     I, expand.grid, unname
#> 
#> 
#> Attaching package: 'Biostrings'
#> 
#> The following object is masked from 'package:base':
#> 
#>     strsplit
#> 
#> 
#> Attaching package: 'rtracklayer'
#> 
#> The following object is masked from 'package:BiocIO':
#> 
#>     FileForFormat
#> 
#> ✔ [2024-09-23 09:59:26.947172]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:26.950115]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:26.951423]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:26.953608]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:26.95429]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:26.956752]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:26.958477]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:26.959509]: SBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:26.96043]: Extracting 5' and 3' adjacent bases.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:27.449461]: Extracting +/- 20bp around mutated bases for background C>T estimation.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:27.661202]: Estimating APOBEC enrichment scores.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:27.662117]: Performing one-way Fisher's test for APOBEC enrichment.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.665693]: APOBEC related mutations are enriched in 0% of samples (APOBEC enrichment score > 2; 0 of 1 samples)
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:27.666988]: Creating SBS sample-by-component matrices.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.66881]: SBS-6 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.670543]: SBS-96 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.675757]: SBS-1536 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:27.676423]: Return SBS-96 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.677512]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:27.678409]: 3.945 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:27.679086]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:27.683314]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.684465]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.68692]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.688253]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.690385]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.691055]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.6935]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.695684]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:27.696564]: DBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:27.697347]: Searching DBS records...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.708653]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.849969]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.855844]: DBS-78 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.861426]: DBS-1248 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:27.862314]: Return SBS-78 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.863309]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:27.864063]: 0.185 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:27.864707]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:27.869171]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.870321]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.876345]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.884074]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.892931]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.896654]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:27.903395]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.292134]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.293832]: INDEL matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.513155]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.514224]: INDEL length extracted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.51528]: Adjacent copies counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.51683]: Microhomology size calculated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.520966]: INDEL records classified into different components (types).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.523018]: ID-28 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.524825]: ID-83 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.525434]: Return ID-83 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.526095]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.526712]: 0.662 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Tally (CopyNumber) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.528668]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.529266]: Genome build  : hg19.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.529918]: Genome measure: wg.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.530506]: When add_loh is TRUE, use_all is forced to TRUE.
#> Please drop columns you don't want to keep before reading.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'sigminer version 2.3.2
#> - Star me at https://github.com/ShixiangWang/sigminer
#> - Run hello() to see usage and citation.
#> ✔ [2024-09-23 09:59:28.535927]: Chromosome size database for build obtained.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.536649]: Reading input.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.537259]: A data frame as input detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.537897]: Column names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.538624]: Column order set.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.539281]: Rows with NA copy number removed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.540523]: Chromosomes unified.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.541841]: Value 2 (normal copy) filled to uncalled chromosomes.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.543257]: Data imported.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.543854]: Segments info:
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.544514]:     Keep - 253
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.545086]:     Drop - 0
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.545938]: Segments sorted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.546531]: Adding LOH labels...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.547303]: Skipped joining adjacent segments with same copy number value.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.547929]: Segmental table cleaned.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.54854]: Annotating.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.55435]: Annotation done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.555005]: Summarizing per sample.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.561313]: Summarized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.561992]: Generating CopyNumber object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.562945]: Generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.563579]: Validating object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.564214]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.564837]: 0.036 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.565457]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.566811]: When you use method 'S', please make sure you have set 'join_adj_seg' to FALSE and 'add_loh' to TRUE in 'read_copynumber() in the previous step!
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.574175]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.57492]: 0.009 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.575573]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.576521]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> Attaching package: 'purrr'
#> 
#> The following object is masked from 'package:XVector':
#> 
#>     compact
#> 
#> The following object is masked from 'package:GenomicRanges':
#> 
#>     reduce
#> 
#> The following object is masked from 'package:IRanges':
#> 
#>     reduce
#> 
#> ℹ [2024-09-23 09:59:28.61032]: Getting breakpoint count per 10 Mb...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.642691]: Getting breakpoint count per chromosome arm...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.670892]: Getting copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.673443]: Getting change-point copy number change...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.680558]: Getting length of chains of oscillating copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.687803]: Getting (log10 based) segment size...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.689647]: Getting the minimal number of chromosome with 50% CNV...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.694536]: Getting burden of chromosome...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.709863]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.710699]: Step: generating copy number components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.711278]: `feature_setting` checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.712575]: Step: counting components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.807581]: Counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.808655]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.810165]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.810939]: 0.235 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.811568]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.812184]: Generated 'Index' column to track the copy number segment location.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.813109]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.815499]: NOTE: this method derives features for each segment. Be patient...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.823818]: Getting absolute copy number value of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.832112]: Getting segment size of eash segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.834337]: Getting context change shape based on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.843672]: Getting change extent on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.852542]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.853516]: Step: generating copy number components based on combination.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.855005]: Classified and combined.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.855634]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.873057]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.874026]: 0.062 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Fitting ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Single Base Substitutions (SBS96) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.877631]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.878289]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.878919]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.879528]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.881478]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.882258]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.882937]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.883731]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.884662]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.885337]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.886051]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.886763]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.887457]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.888152]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.904149]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.905432]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.906378]: Fitting sample: DO1000
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.909703]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.910637]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.916285]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.917086]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:28.917874]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.918553]: 0.037 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.92022]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.920888]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:28.931412]: Processing sample `DO1000`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.6552]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.656208]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.656879]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.657425]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.662299]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.662921]: 0.006 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.663552]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.664097]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.669194]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.67008]: Total 1.792 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── INDELS (ID83) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.671377]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.671965]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.672503]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.673036]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.674298]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.674842]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.675385]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.675976]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.676647]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.677192]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.677752]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.678281]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.678813]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.679344]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.679931]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.680466]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.681072]: Fitting sample: DO1000
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.681706]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.682238]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.685532]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.686107]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:30.686709]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.687265]: 0.013 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.688401]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:30.688948]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:32.108293]: Processing sample `DO1000`.
#> ✔ [2024-09-23 09:59:33.383465]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.384327]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.38497]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:33.385523]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:33.386934]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.387498]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:33.38805]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.388585]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:33.392845]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.393621]: Total 2.722 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Doublet Mutations (DBS78) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.395347]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.397519]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.400366]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.401793]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.403327]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.403928]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:33.404563]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:33.405256]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:33.405989]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.406617]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:33.407881]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.408646]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:33.409316]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:33.409952]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:33.410586]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.411217]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.412161]: Fitting sample: DO1000
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:33.412967]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.413604]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:33.420063]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.420968]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:33.421623]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.422198]: 0.019 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.423393]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:33.423941]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:34.356345]: Processing sample `DO1000`.
#> ✔ [2024-09-23 09:59:35.690674]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.691718]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.692351]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:35.692905]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:35.694215]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.694796]: 0.002 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:35.695345]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.695872]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:35.700347]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.701178]: Total 2.306 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Copy Number Alterations (CN48) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.703452]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.704233]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.705638]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.706805]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.709683]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.710393]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:35.711076]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:35.711809]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:35.712562]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.713415]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:35.71415]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.714818]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:35.715456]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:35.716168]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:35.717757]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.718807]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.720108]: Fitting sample: DO1000
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:35.721587]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.722446]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:35.732265]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.733367]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:35.734425]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.73515]: 0.025 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.737013]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:35.737925]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:36.80707]: Processing sample `DO1000`.
#> ✔ [2024-09-23 09:59:38.175787]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:38.176777]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:38.177448]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:38.178045]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:38.179703]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:38.180314]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:38.180904]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:38.181469]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-23 09:59:38.186081]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-23 09:59:38.186919]: Total 2.483 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Write Output ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Raw counts (tally) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ SBS96 tally written to csv: 'pcawg_signature_results/DO1000/SBS96_catalogue.DO1000.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ SBS1536 tally written to csv: 'pcawg_signature_results/DO1000/SBS1536_catalogue.DO1000.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ ID83 tally written to csv: 'pcawg_signature_results/DO1000/ID83_catalogue.DO1000.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ DBS78 tally written to csv: 'pcawg_signature_results/DO1000/DBS78_catalogue.DO1000.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ CN48 tally written to csv: 'pcawg_signature_results/DO1000/CN48_catalogue.DO1000.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ CN80 tally written to csv: 'pcawg_signature_results/DO1000/CN80_catalogue.DO1000.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ CN176 tally written to csv: 'pcawg_signature_results/DO1000/CN176_catalogue.DO1000.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Fit (Exposures) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ SBS96 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1000/SBS96_fit.DO1000.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ SBS96 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1000/SBS96_fit.DO1000.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ SBS96 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1000/SBS96_fit.DO1000.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ SBS96 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1000/SBS96_fit.DO1000.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ SBS96 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1000/SBS96_fit.DO1000.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ SBS96 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1000/SBS96_fit.DO1000.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ DBS78 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1000/DBS78_fit.DO1000.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ DBS78 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1000/DBS78_fit.DO1000.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ DBS78 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1000/DBS78_fit.DO1000.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ DBS78 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1000/DBS78_fit.DO1000.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ DBS78 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1000/DBS78_fit.DO1000.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ DBS78 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1000/DBS78_fit.DO1000.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ ID83 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1000/ID83_fit.DO1000.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ ID83 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1000/ID83_fit.DO1000.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ ID83 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1000/ID83_fit.DO1000.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ ID83 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1000/ID83_fit.DO1000.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ ID83 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1000/ID83_fit.DO1000.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ ID83 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1000/ID83_fit.DO1000.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ CN48 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1000/CN48_fit.DO1000.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ CN48 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1000/CN48_fit.DO1000.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ CN48 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1000/CN48_fit.DO1000.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ CN48 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1000/CN48_fit.DO1000.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ CN48 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1000/CN48_fit.DO1000.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ CN48 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1000/CN48_fit.DO1000.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#>  Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100%
#> ✔ Creating Output Directory at 'pcawg_signature_results/DO1000' [16.5s]
#> 
#> ℹ Finished successfully
#> ✔ Finished successfully [11ms]
#> 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Found 1 samples described in the VCF [DO1001]
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ Returning data for only sample [DO1001] samples in format column of VCF.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Running Signature Analysis. This will take some time ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Mutational Signature Analysis ───────────────────────────────────────────────
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Checking arguments ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ Output Directory: 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ Reference Genome: BSgenome.Hsapiens.UCSC.hg19
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'Please note that the generated MAF object is designed for mutational signature analysis, not recommended for Maftools analysis!
#> 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Tally (Small Variants) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:39.332446]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:39.336902]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:39.337671]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:39.339947]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:39.342621]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:39.34899]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:39.349728]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:39.353894]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:39.355456]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:39.356013]: SBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:39.356843]: Extracting 5' and 3' adjacent bases.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.051134]: Extracting +/- 20bp around mutated bases for background C>T estimation.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.269323]: Estimating APOBEC enrichment scores.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.270244]: Performing one-way Fisher's test for APOBEC enrichment.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.272372]: APOBEC related mutations are enriched in 0% of samples (APOBEC enrichment score > 2; 0 of 1 samples)
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.27357]: Creating SBS sample-by-component matrices.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.275258]: SBS-6 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.277008]: SBS-96 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.282413]: SBS-1536 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.283062]: Return SBS-96 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.284136]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.284757]: 0.952 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.285383]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.288782]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.289489]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.291816]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.294265]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.299612]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.306635]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.311474]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.313165]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.313726]: DBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.314464]: Searching DBS records...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.335309]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.487733]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.491799]: DBS-78 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.496051]: DBS-1248 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.496691]: Return SBS-78 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.497575]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.498155]: 0.213 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.498743]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.501929]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.502643]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.504944]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.507334]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.512664]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.513441]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.518367]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.520174]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.520758]: INDEL matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.800377]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.801385]: INDEL length extracted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.802916]: Adjacent copies counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.814182]: Microhomology size calculated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.820341]: INDEL records classified into different components (types).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.822221]: ID-28 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.823841]: ID-83 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.82441]: Return ID-83 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.825015]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.825585]: 0.327 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Tally (CopyNumber) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.827003]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.827557]: Genome build  : hg19.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.828103]: Genome measure: wg.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.828644]: When add_loh is TRUE, use_all is forced to TRUE.
#> Please drop columns you don't want to keep before reading.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.830164]: Chromosome size database for build obtained.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.830726]: Reading input.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.831284]: A data frame as input detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.831876]: Column names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.832456]: Column order set.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.833067]: Rows with NA copy number removed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.834583]: Chromosomes unified.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.835775]: Value 2 (normal copy) filled to uncalled chromosomes.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.837093]: Data imported.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.837648]: Segments info:
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.838185]:     Keep - 637
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.838734]:     Drop - 0
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.839414]: Segments sorted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.839975]: Adding LOH labels...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.840714]: Skipped joining adjacent segments with same copy number value.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.841318]: Segmental table cleaned.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.841882]: Annotating.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.847442]: Annotation done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.848047]: Summarizing per sample.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.852413]: Summarized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.853028]: Generating CopyNumber object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.853691]: Generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.854251]: Validating object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.854799]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.855366]: 0.028 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.855931]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.856795]: When you use method 'S', please make sure you have set 'join_adj_seg' to FALSE and 'add_loh' to TRUE in 'read_copynumber() in the previous step!
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.861501]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.862115]: 0.006 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.862685]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.863524]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.875578]: Getting breakpoint count per 10 Mb...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.889508]: Getting breakpoint count per chromosome arm...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.896268]: Getting copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.897799]: Getting change-point copy number change...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.903553]: Getting length of chains of oscillating copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.910105]: Getting (log10 based) segment size...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.911641]: Getting the minimal number of chromosome with 50% CNV...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.91533]: Getting burden of chromosome...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.918657]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.919271]: Step: generating copy number components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:40.919851]: `feature_setting` checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:40.921144]: Step: counting components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.01225]: Counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.01309]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.0139]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.01449]: 0.152 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.015059]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.015643]: Generated 'Index' column to track the copy number segment location.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.016522]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.018559]: NOTE: this method derives features for each segment. Be patient...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.026046]: Getting absolute copy number value of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.032613]: Getting segment size of eash segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.034168]: Getting context change shape based on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.041146]: Getting change extent on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.050774]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.051453]: Step: generating copy number components based on combination.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.052458]: Classified and combined.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.053024]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.064874]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.06557]: 0.05 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Fitting ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Single Base Substitutions (SBS96) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.068498]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.069055]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.069641]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.070183]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.071453]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.071997]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.072548]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.073241]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.074064]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.074603]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.075141]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.075681]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.076237]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.076774]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.077314]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.07787]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.078497]: Fitting sample: DO1001
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.080171]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.080716]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.084002]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.084569]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:41.085164]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.085729]: 0.014 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.086834]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.087383]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:41.09866]: Processing sample `DO1001`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.846101]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.846886]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.847506]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.848057]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.851463]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.85205]: 0.005 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.852643]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.853194]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.858182]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.859087]: Total 1.791 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── INDELS (ID83) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.860427]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.86099]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.86154]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.86208]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.863343]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.863884]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.864437]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.865041]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.865683]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.866232]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.866771]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.867316]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.867866]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.868406]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.868971]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.869508]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.870125]: Fitting sample: DO1001
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.870757]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.871293]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.874479]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.875047]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:42.875632]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.876189]: 0.013 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.877278]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:42.877829]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:43.808529]: Processing sample `DO1001`.
#> ✔ [2024-09-23 09:59:45.098258]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.099091]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.099722]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:45.100276]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:45.101689]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.10227]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:45.102829]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.103365]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:45.107688]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.108508]: Total 2.248 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Doublet Mutations (DBS78) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.110666]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.111482]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.112355]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.114885]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.116703]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.117736]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:45.11844]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:45.119164]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:45.119931]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.12059]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:45.121251]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.121898]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:45.123745]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:45.124508]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:45.125543]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.12645]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.1273]: Fitting sample: DO1001
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:45.128094]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.128771]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:45.133636]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.13447]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:45.135092]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.135662]: 0.019 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.136836]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:45.137401]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:46.168684]: Processing sample `DO1001`.
#> ✔ [2024-09-23 09:59:47.602316]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.603239]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.603958]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:47.604629]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:47.606094]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.606718]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:47.607345]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.607935]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:47.612528]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.613378]: Total 2.503 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Copy Number Alterations (CN48) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.615687]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.618897]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.621905]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.622765]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.624407]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.62509]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:47.625786]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:47.626528]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:47.627319]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.628107]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:47.630393]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.631228]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:47.631967]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:47.63267]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:47.633358]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.634034]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.634856]: Fitting sample: DO1001
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:47.635804]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.637105]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:47.642594]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.643399]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:47.644051]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.64473]: 0.02 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.645987]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:47.646585]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:48.616942]: Processing sample `DO1001`.
#> ✔ [2024-09-23 09:59:50.134896]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:50.135826]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:50.136556]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:50.143836]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:50.145438]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:50.146033]: 0.009 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:50.146704]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:50.147287]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-23 09:59:50.151724]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-23 09:59:50.15258]: Total 2.537 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Write Output ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Raw counts (tally) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ SBS96 tally written to csv: 'pcawg_signature_results/DO1001/SBS96_catalogue.DO1001.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ SBS1536 tally written to csv: 'pcawg_signature_results/DO1001/SBS1536_catalogue.DO1001.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ ID83 tally written to csv: 'pcawg_signature_results/DO1001/ID83_catalogue.DO1001.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ DBS78 tally written to csv: 'pcawg_signature_results/DO1001/DBS78_catalogue.DO1001.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ CN48 tally written to csv: 'pcawg_signature_results/DO1001/CN48_catalogue.DO1001.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ CN80 tally written to csv: 'pcawg_signature_results/DO1001/CN80_catalogue.DO1001.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ CN176 tally written to csv: 'pcawg_signature_results/DO1001/CN176_catalogue.DO1001.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Fit (Exposures) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ SBS96 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1001/SBS96_fit.DO1001.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ SBS96 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1001/SBS96_fit.DO1001.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ SBS96 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1001/SBS96_fit.DO1001.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ SBS96 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1001/SBS96_fit.DO1001.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ SBS96 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1001/SBS96_fit.DO1001.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ SBS96 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1001/SBS96_fit.DO1001.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ DBS78 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1001/DBS78_fit.DO1001.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ DBS78 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1001/DBS78_fit.DO1001.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ DBS78 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1001/DBS78_fit.DO1001.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ DBS78 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1001/DBS78_fit.DO1001.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ DBS78 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1001/DBS78_fit.DO1001.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ DBS78 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1001/DBS78_fit.DO1001.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ ID83 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1001/ID83_fit.DO1001.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ ID83 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1001/ID83_fit.DO1001.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ ID83 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1001/ID83_fit.DO1001.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ ID83 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1001/ID83_fit.DO1001.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ ID83 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1001/ID83_fit.DO1001.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ ID83 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1001/ID83_fit.DO1001.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ CN48 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1001/CN48_fit.DO1001.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ CN48 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1001/CN48_fit.DO1001.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ CN48 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1001/CN48_fit.DO1001.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ CN48 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1001/CN48_fit.DO1001.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ CN48 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1001/CN48_fit.DO1001.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ CN48 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1001/CN48_fit.DO1001.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#>  Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100%
#> ✔ Creating Output Directory at 'pcawg_signature_results/DO1001' [12.5s]
#> 
#> ℹ Finished successfully
#> ✔ Finished successfully [7ms]
#> 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Found 1 samples described in the VCF [DO1002]
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ Returning data for only sample [DO1002] samples in format column of VCF.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Running Signature Analysis. This will take some time ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Mutational Signature Analysis ───────────────────────────────────────────────
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Checking arguments ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ Output Directory: 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ Reference Genome: BSgenome.Hsapiens.UCSC.hg19
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'Please note that the generated MAF object is designed for mutational signature analysis, not recommended for Maftools analysis!
#> 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Tally (Small Variants) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:51.819418]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:51.822907]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:51.823664]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:51.825968]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:51.82838]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:51.833788]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:51.83456]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:51.838948]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:51.840506]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:51.841068]: SBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:51.841913]: Extracting 5' and 3' adjacent bases.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.450615]: Extracting +/- 20bp around mutated bases for background C>T estimation.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.658693]: Estimating APOBEC enrichment scores.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.659616]: Performing one-way Fisher's test for APOBEC enrichment.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.661728]: APOBEC related mutations are enriched in 0% of samples (APOBEC enrichment score > 2; 0 of 1 samples)
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.662937]: Creating SBS sample-by-component matrices.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.664603]: SBS-6 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.666336]: SBS-96 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.671458]: SBS-1536 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.672057]: Return SBS-96 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.673041]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.673627]: 0.854 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.67422]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.677314]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.678011]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.680294]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.682738]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.688118]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.6889]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.693444]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.695051]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.695619]: DBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.696376]: Searching DBS records...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.716964]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.943606]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.947968]: DBS-78 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.95792]: DBS-1248 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.958659]: Return SBS-78 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.959568]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.960168]: 0.286 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.960809]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.963886]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.964582]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.966858]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.969632]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.976738]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.97751]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.982047]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:52.983826]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:52.984395]: INDEL matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.286494]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.287644]: INDEL length extracted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.289793]: Adjacent copies counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.314684]: Microhomology size calculated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.323899]: INDEL records classified into different components (types).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.325679]: ID-28 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.327362]: ID-83 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.327949]: Return ID-83 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.32858]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.329158]: 0.368 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Tally (CopyNumber) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.330611]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.331174]: Genome build  : hg19.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.331724]: Genome measure: wg.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.332274]: When add_loh is TRUE, use_all is forced to TRUE.
#> Please drop columns you don't want to keep before reading.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.333786]: Chromosome size database for build obtained.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.334384]: Reading input.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.334958]: A data frame as input detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.335555]: Column names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.336147]: Column order set.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.336781]: Rows with NA copy number removed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.338251]: Chromosomes unified.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.339436]: Value 2 (normal copy) filled to uncalled chromosomes.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.340758]: Data imported.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.341325]: Segments info:
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.341879]:     Keep - 573
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.34246]:     Drop - 0
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.343136]: Segments sorted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.343683]: Adding LOH labels...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.344422]: Skipped joining adjacent segments with same copy number value.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.34502]: Segmental table cleaned.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.345565]: Annotating.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.350928]: Annotation done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.351529]: Summarizing per sample.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.355833]: Summarized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.356447]: Generating CopyNumber object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.357155]: Generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.357716]: Validating object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.358272]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.358871]: 0.028 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.359433]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.360314]: When you use method 'S', please make sure you have set 'join_adj_seg' to FALSE and 'add_loh' to TRUE in 'read_copynumber() in the previous step!
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.365132]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.365748]: 0.006 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.366326]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.3672]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.379206]: Getting breakpoint count per 10 Mb...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.393004]: Getting breakpoint count per chromosome arm...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.399722]: Getting copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.401219]: Getting change-point copy number change...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.406872]: Getting length of chains of oscillating copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.413148]: Getting (log10 based) segment size...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.41466]: Getting the minimal number of chromosome with 50% CNV...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.418329]: Getting burden of chromosome...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.421724]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.422341]: Step: generating copy number components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.422916]: `feature_setting` checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.424238]: Step: counting components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.508719]: Counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.509601]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.510436]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.511032]: 0.145 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.511625]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.512232]: Generated 'Index' column to track the copy number segment location.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.513155]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.515186]: NOTE: this method derives features for each segment. Be patient...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.523192]: Getting absolute copy number value of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.529845]: Getting segment size of eash segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.531407]: Getting context change shape based on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.538402]: Getting change extent on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.545494]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.546143]: Step: generating copy number components based on combination.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.547143]: Classified and combined.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.547746]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.56012]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.560871]: 0.049 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Fitting ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Single Base Substitutions (SBS96) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.564]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.564586]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.565164]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.565756]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.567066]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.567659]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.568241]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.568965]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.56987]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.570471]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.57105]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.571632]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.5722]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.572765]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.573343]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.573916]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.574577]: Fitting sample: DO1002
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.576393]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.576967]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.58045]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.581067]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:53.581702]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.58231]: 0.015 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.583456]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.584042]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:53.594705]: Processing sample `DO1002`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.487559]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.488612]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.489523]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.490418]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.49441]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.495342]: 0.006 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.495983]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.496652]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.502764]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.503948]: Total 1.94 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── INDELS (ID83) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.506322]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.50699]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.507605]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.508195]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.509569]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.510145]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.510755]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.511503]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.512203]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.512793]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.513368]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.513936]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.51451]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.515081]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.515691]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.516265]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.516974]: Fitting sample: DO1002
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.517676]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.518353]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.525669]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.52709]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:55.528214]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.529128]: 0.019 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.530577]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:55.531247]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:56.659872]: Processing sample `DO1002`.
#> ✔ [2024-09-23 09:59:58.237556]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.238489]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.239154]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:58.239748]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:58.241326]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.241967]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:58.242562]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.243123]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:58.247557]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.248478]: Total 2.742 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Doublet Mutations (DBS78) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.250989]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.255297]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.257608]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.258466]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.260016]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.260669]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:58.26136]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:58.2621]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:58.262883]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.26378]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:58.26504]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.266215]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:58.266953]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:58.269808]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:58.272184]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.274782]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.275675]: Fitting sample: DO1002
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:58.27669]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.27753]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:58.281755]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.282642]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 09:59:58.283374]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.284]: 0.024 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.285179]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:58.285768]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 09:59:59.385361]: Processing sample `DO1002`.
#> ✔ [2024-09-23 10:00:00.802968]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.80401]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.804696]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:00.805301]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:00.806713]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.807386]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:00.807997]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.808582]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:00.81307]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.813965]: Total 2.563 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Copy Number Alterations (CN48) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.815812]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.817226]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.820256]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.822436]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.824185]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.824936]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:00.825676]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:00.826488]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:00.827723]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.828604]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:00.829363]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.830083]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:00.830794]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:00.831497]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:00.832364]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.833407]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.834341]: Fitting sample: DO1002
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:00.835905]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.837231]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:00.842572]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.843417]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:00.84409]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.844722]: 0.021 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.845966]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:00.846578]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:01.869485]: Processing sample `DO1002`.
#> ✔ [2024-09-23 10:00:03.217151]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:03.218039]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:03.218695]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:03.219305]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:03.220956]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:03.221604]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:03.222227]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:03.222818]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-23 10:00:03.227325]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-23 10:00:03.228268]: Total 2.412 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Write Output ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Raw counts (tally) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ SBS96 tally written to csv: 'pcawg_signature_results/DO1002/SBS96_catalogue.DO1002.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ SBS1536 tally written to csv: 'pcawg_signature_results/DO1002/SBS1536_catalogue.DO1002.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ ID83 tally written to csv: 'pcawg_signature_results/DO1002/ID83_catalogue.DO1002.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ DBS78 tally written to csv: 'pcawg_signature_results/DO1002/DBS78_catalogue.DO1002.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ CN48 tally written to csv: 'pcawg_signature_results/DO1002/CN48_catalogue.DO1002.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ CN80 tally written to csv: 'pcawg_signature_results/DO1002/CN80_catalogue.DO1002.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ CN176 tally written to csv: 'pcawg_signature_results/DO1002/CN176_catalogue.DO1002.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Fit (Exposures) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ SBS96 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1002/SBS96_fit.DO1002.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ SBS96 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1002/SBS96_fit.DO1002.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ SBS96 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1002/SBS96_fit.DO1002.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ SBS96 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1002/SBS96_fit.DO1002.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ SBS96 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1002/SBS96_fit.DO1002.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ SBS96 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1002/SBS96_fit.DO1002.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ DBS78 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1002/DBS78_fit.DO1002.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ DBS78 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1002/DBS78_fit.DO1002.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ DBS78 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1002/DBS78_fit.DO1002.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ DBS78 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1002/DBS78_fit.DO1002.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ DBS78 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1002/DBS78_fit.DO1002.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ DBS78 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1002/DBS78_fit.DO1002.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ ID83 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1002/ID83_fit.DO1002.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ ID83 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1002/ID83_fit.DO1002.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ ID83 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1002/ID83_fit.DO1002.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ ID83 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1002/ID83_fit.DO1002.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ ID83 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1002/ID83_fit.DO1002.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ ID83 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1002/ID83_fit.DO1002.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ CN48 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1002/CN48_fit.DO1002.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ CN48 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1002/CN48_fit.DO1002.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ CN48 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1002/CN48_fit.DO1002.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ CN48 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1002/CN48_fit.DO1002.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ CN48 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1002/CN48_fit.DO1002.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ CN48 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1002/CN48_fit.DO1002.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#>  Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100%
#> ✔ Creating Output Directory at 'pcawg_signature_results/DO1002' [12.5s]
#> 
#> ℹ Finished successfully
#> ✔ Finished successfully [7ms]
#> 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Found 1 samples described in the VCF [DO1003]
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ Returning data for only sample [DO1003] samples in format column of VCF.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Running Signature Analysis. This will take some time ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Mutational Signature Analysis ───────────────────────────────────────────────
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Checking arguments ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ Output Directory: 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ Reference Genome: BSgenome.Hsapiens.UCSC.hg19
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'Please note that the generated MAF object is designed for mutational signature analysis, not recommended for Maftools analysis!
#> 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Tally (Small Variants) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:04.30161]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:04.306099]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:04.306954]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:04.309202]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:04.311495]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:04.316269]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:04.31707]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:04.321121]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:04.322935]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:04.323581]: SBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:04.324475]: Extracting 5' and 3' adjacent bases.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:04.908018]: Extracting +/- 20bp around mutated bases for background C>T estimation.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.067078]: Estimating APOBEC enrichment scores.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.067998]: Performing one-way Fisher's test for APOBEC enrichment.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.069921]: APOBEC related mutations are enriched in 0% of samples (APOBEC enrichment score > 2; 0 of 1 samples)
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.071149]: Creating SBS sample-by-component matrices.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.072718]: SBS-6 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.074432]: SBS-96 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.079379]: SBS-1536 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.080004]: Return SBS-96 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.081029]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.081623]: 0.78 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.082215]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.085256]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.085955]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.088094]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.090244]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.094796]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.095537]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.099307]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.100856]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.101445]: DBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.102209]: Searching DBS records...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.122626]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.331857]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.337383]: DBS-78 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.342353]: DBS-1248 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.343839]: Return SBS-78 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.345148]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.34595]: 0.264 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.346667]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.349982]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.350802]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.35369]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.356143]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.361224]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.363015]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.367653]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.369828]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.370483]: INDEL matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.720706]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.721782]: INDEL length extracted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.723065]: Adjacent copies counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.73164]: Microhomology size calculated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.737324]: INDEL records classified into different components (types).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.744586]: ID-28 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.752978]: ID-83 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.755274]: Return ID-83 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.758538]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.760941]: 0.414 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Tally (CopyNumber) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.767156]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.768576]: Genome build  : hg19.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.770698]: Genome measure: wg.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.772757]: When add_loh is TRUE, use_all is forced to TRUE.
#> Please drop columns you don't want to keep before reading.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.777231]: Chromosome size database for build obtained.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.78039]: Reading input.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.784122]: A data frame as input detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.7867]: Column names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.787687]: Column order set.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.788504]: Rows with NA copy number removed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.790812]: Chromosomes unified.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.792792]: Value 2 (normal copy) filled to uncalled chromosomes.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.795009]: Data imported.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.795936]: Segments info:
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.796676]:     Keep - 1007
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.798058]:     Drop - 0
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.799591]: Segments sorted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.801916]: Adding LOH labels...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.805056]: Skipped joining adjacent segments with same copy number value.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.807102]: Segmental table cleaned.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.8081]: Annotating.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.815733]: Annotation done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.816505]: Summarizing per sample.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.821245]: Summarized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.821951]: Generating CopyNumber object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.822739]: Generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.823345]: Validating object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.823954]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.824559]: 0.057 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.825175]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.826104]: When you use method 'S', please make sure you have set 'join_adj_seg' to FALSE and 'add_loh' to TRUE in 'read_copynumber() in the previous step!
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.831341]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.832002]: 0.007 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.832607]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.833481]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.845845]: Getting breakpoint count per 10 Mb...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.866975]: Getting breakpoint count per chromosome arm...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.874468]: Getting copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.876146]: Getting change-point copy number change...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.882249]: Getting length of chains of oscillating copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.889365]: Getting (log10 based) segment size...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.89098]: Getting the minimal number of chromosome with 50% CNV...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.894751]: Getting burden of chromosome...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.898563]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.899207]: Step: generating copy number components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.899799]: `feature_setting` checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.90109]: Step: counting components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.988638]: Counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.989638]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.990479]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.991098]: 0.158 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.991709]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:05.992334]: Generated 'Index' column to track the copy number segment location.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.99328]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:05.995371]: NOTE: this method derives features for each segment. Be patient...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.003126]: Getting absolute copy number value of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.009826]: Getting segment size of eash segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.011416]: Getting context change shape based on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.01879]: Getting change extent on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:06.02671]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.027432]: Step: generating copy number components based on combination.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:06.028512]: Classified and combined.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.029182]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:06.043447]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.04434]: 0.053 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Fitting ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Single Base Substitutions (SBS96) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.047649]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.048302]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.048933]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.049589]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.051015]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.051634]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:06.05226]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:06.053056]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:06.053961]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.054586]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:06.055196]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.055849]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:06.056464]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:06.057105]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:06.057722]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.058381]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.05909]: Fitting sample: DO1003
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:06.061011]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.061591]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:06.065291]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.066021]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:06.066719]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.06735]: 0.016 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.068602]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.069237]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:06.079814]: Processing sample `DO1003`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.887639]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.888463]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.8891]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.889669]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.893037]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.893665]: 0.005 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.894252]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.894808]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.899777]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.900687]: Total 1.853 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── INDELS (ID83) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.902071]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.902648]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.903221]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.903782]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.905066]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.905627]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.906207]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.906837]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.907499]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.908076]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.90864]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.909199]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.909772]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.910338]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.910902]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.91146]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.912093]: Fitting sample: DO1003
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.912755]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.913317]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.916469]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.917061]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:07.917677]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.91826]: 0.013 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.919341]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:07.919916]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:08.878915]: Processing sample `DO1003`.
#> ✔ [2024-09-23 10:00:10.146966]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.148006]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.148726]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:10.149355]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:10.150838]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.151481]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:10.152126]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.152724]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:10.157358]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.158324]: Total 2.256 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Doublet Mutations (DBS78) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.162611]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.163883]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.164818]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.16565]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.167944]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.168761]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:10.169513]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:10.170311]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:10.171585]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.17293]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:10.173814]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.174562]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:10.175381]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:10.176764]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:10.177638]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.178403]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.179246]: Fitting sample: DO1003
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:10.180271]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.181077]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:10.187364]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.188367]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:10.189187]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.18993]: 0.022 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.191415]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:10.19216]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:11.13563]: Processing sample `DO1003`.
#> ✔ [2024-09-23 10:00:12.419135]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.420051]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.420785]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:12.421423]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:12.422908]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.423577]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:12.424206]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.424805]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:12.42934]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.430203]: Total 2.268 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Copy Number Alterations (CN48) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.432129]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.432868]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.433576]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.434298]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.436052]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.437485]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:12.438426]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:12.439307]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:12.440188]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.440941]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:12.441688]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.442434]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:12.443173]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:12.443891]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:12.444588]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.445589]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.446655]: Fitting sample: DO1003
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:12.447747]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.448609]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:12.453684]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.455237]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:12.456216]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.457086]: 0.021 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.458503]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:12.459151]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:13.400512]: Processing sample `DO1003`.
#> ✔ [2024-09-23 10:00:14.679101]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:14.680032]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:14.680735]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:14.681384]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:14.683208]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:14.683949]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:14.684634]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:14.685298]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-23 10:00:14.690156]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-23 10:00:14.691188]: Total 2.259 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Write Output ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Raw counts (tally) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ SBS96 tally written to csv: 'pcawg_signature_results/DO1003/SBS96_catalogue.DO1003.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ SBS1536 tally written to csv: 'pcawg_signature_results/DO1003/SBS1536_catalogue.DO1003.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ ID83 tally written to csv: 'pcawg_signature_results/DO1003/ID83_catalogue.DO1003.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ DBS78 tally written to csv: 'pcawg_signature_results/DO1003/DBS78_catalogue.DO1003.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ CN48 tally written to csv: 'pcawg_signature_results/DO1003/CN48_catalogue.DO1003.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ CN80 tally written to csv: 'pcawg_signature_results/DO1003/CN80_catalogue.DO1003.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ CN176 tally written to csv: 'pcawg_signature_results/DO1003/CN176_catalogue.DO1003.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Fit (Exposures) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ SBS96 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1003/SBS96_fit.DO1003.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ SBS96 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1003/SBS96_fit.DO1003.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ SBS96 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1003/SBS96_fit.DO1003.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ SBS96 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1003/SBS96_fit.DO1003.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ SBS96 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1003/SBS96_fit.DO1003.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ SBS96 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1003/SBS96_fit.DO1003.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ DBS78 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1003/DBS78_fit.DO1003.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ DBS78 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1003/DBS78_fit.DO1003.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ DBS78 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1003/DBS78_fit.DO1003.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ DBS78 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1003/DBS78_fit.DO1003.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ DBS78 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1003/DBS78_fit.DO1003.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ DBS78 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1003/DBS78_fit.DO1003.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ ID83 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1003/ID83_fit.DO1003.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ ID83 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1003/ID83_fit.DO1003.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ ID83 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1003/ID83_fit.DO1003.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ ID83 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1003/ID83_fit.DO1003.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ ID83 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1003/ID83_fit.DO1003.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ ID83 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1003/ID83_fit.DO1003.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ CN48 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1003/CN48_fit.DO1003.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ CN48 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1003/CN48_fit.DO1003.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ CN48 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1003/CN48_fit.DO1003.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ CN48 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1003/CN48_fit.DO1003.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ CN48 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1003/CN48_fit.DO1003.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ CN48 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1003/CN48_fit.DO1003.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#>  Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100%
#> ✔ Creating Output Directory at 'pcawg_signature_results/DO1003' [11.5s]
#> 
#> ℹ Finished successfully
#> ✔ Finished successfully [6ms]
#> 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Found 1 samples described in the VCF [DO1004]
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ Returning data for only sample [DO1004] samples in format column of VCF.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Running Signature Analysis. This will take some time ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Mutational Signature Analysis ───────────────────────────────────────────────
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Checking arguments ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ Output Directory: 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ Reference Genome: BSgenome.Hsapiens.UCSC.hg19
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'Please note that the generated MAF object is designed for mutational signature analysis, not recommended for Maftools analysis!
#> 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Tally (Small Variants) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:15.934128]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:15.937895]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:15.938675]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:15.941826]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:15.946171]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:15.956579]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:15.957542]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:15.96599]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:15.967754]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:15.968357]: SBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:15.969382]: Extracting 5' and 3' adjacent bases.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:16.730674]: Extracting +/- 20bp around mutated bases for background C>T estimation.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:16.918806]: Estimating APOBEC enrichment scores.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:16.9198]: Performing one-way Fisher's test for APOBEC enrichment.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:16.928743]: APOBEC related mutations are enriched in 0% of samples (APOBEC enrichment score > 2; 0 of 1 samples)
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:16.930125]: Creating SBS sample-by-component matrices.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:16.931901]: SBS-6 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:16.933761]: SBS-96 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:16.939246]: SBS-1536 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:16.939925]: Return SBS-96 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:16.941022]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:16.941674]: 1.008 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:16.942328]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:16.945471]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:16.946237]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:16.949216]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:16.953618]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:16.964711]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:16.965651]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:16.974591]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:16.976635]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:16.977233]: DBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:16.978228]: Searching DBS records...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.018649]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.257009]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.261501]: DBS-78 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.266066]: DBS-1248 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.2668]: Return SBS-78 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.267811]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.268482]: 0.326 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.269152]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.272451]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.273219]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.276506]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.280901]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.298678]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.299927]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.3089]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.310747]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.311358]: INDEL matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.63138]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.632519]: INDEL length extracted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.635177]: Adjacent copies counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.666176]: Microhomology size calculated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.678006]: INDEL records classified into different components (types).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.679835]: ID-28 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.681521]: ID-83 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.682139]: Return ID-83 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.682797]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.6834]: 0.414 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Tally (CopyNumber) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.684975]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.685574]: Genome build  : hg19.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.68616]: Genome measure: wg.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.686756]: When add_loh is TRUE, use_all is forced to TRUE.
#> Please drop columns you don't want to keep before reading.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.688302]: Chromosome size database for build obtained.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.688906]: Reading input.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.689509]: A data frame as input detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.690132]: Column names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.690764]: Column order set.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.691406]: Rows with NA copy number removed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.692774]: Chromosomes unified.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.694014]: Value 2 (normal copy) filled to uncalled chromosomes.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.695368]: Data imported.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.695979]: Segments info:
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.69656]:     Keep - 442
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.697147]:     Drop - 0
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.697847]: Segments sorted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.698433]: Adding LOH labels...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.699202]: Skipped joining adjacent segments with same copy number value.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.699825]: Segmental table cleaned.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.700409]: Annotating.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.705582]: Annotation done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.706199]: Summarizing per sample.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.710513]: Summarized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.711144]: Generating CopyNumber object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.711832]: Generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.712423]: Validating object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.713004]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.713594]: 0.029 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.714213]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.715095]: When you use method 'S', please make sure you have set 'join_adj_seg' to FALSE and 'add_loh' to TRUE in 'read_copynumber() in the previous step!
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.719855]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.720496]: 0.006 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.721099]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.721964]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.73599]: Getting breakpoint count per 10 Mb...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.749785]: Getting breakpoint count per chromosome arm...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.756443]: Getting copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.757968]: Getting change-point copy number change...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.763733]: Getting length of chains of oscillating copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.769873]: Getting (log10 based) segment size...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.771396]: Getting the minimal number of chromosome with 50% CNV...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.775086]: Getting burden of chromosome...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.778754]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.77939]: Step: generating copy number components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.77999]: `feature_setting` checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.7813]: Step: counting components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.867029]: Counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.867953]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.868831]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.869474]: 0.148 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.870095]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.870761]: Generated 'Index' column to track the copy number segment location.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.87168]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.873801]: NOTE: this method derives features for each segment. Be patient...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.883783]: Getting absolute copy number value of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.890471]: Getting segment size of eash segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.892044]: Getting context change shape based on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.899241]: Getting change extent on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.906533]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.907231]: Step: generating copy number components based on combination.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.908193]: Classified and combined.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.908805]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.920864]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.92159]: 0.051 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Fitting ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Single Base Substitutions (SBS96) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.931351]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.931968]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.932571]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.933154]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.934462]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.93504]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.935633]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.936381]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.937236]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.937816]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.938412]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.939032]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.939621]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.940193]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.940773]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.941347]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.94202]: Fitting sample: DO1004
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.943864]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.944441]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.947975]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.948585]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:17.949217]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.949816]: 0.015 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.950963]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.951548]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:17.963838]: Processing sample `DO1004`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.736339]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.7372]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.737864]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.738451]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.741905]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.74257]: 0.005 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.743159]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.743725]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.748651]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.749564]: Total 1.818 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── INDELS (ID83) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.750977]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.751579]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.752157]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.75273]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.754019]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.754593]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.755177]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.755811]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.756485]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.757076]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.757665]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.758235]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.758805]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.759413]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.760017]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.760589]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.761238]: Fitting sample: DO1004
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.761912]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.762482]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.765777]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.766393]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:19.767016]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.767608]: 0.014 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.768759]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:19.769347]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:20.747934]: Processing sample `DO1004`.
#> ✔ [2024-09-23 10:00:22.047496]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.048437]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.049113]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:22.049713]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:22.051157]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.051784]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:22.052382]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.052958]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:22.057288]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.058147]: Total 2.307 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Doublet Mutations (DBS78) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.060012]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.060733]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.061419]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.062084]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.06384]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.064563]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:22.066941]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:22.068034]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:22.068953]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.069731]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:22.070429]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.071117]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:22.071796]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:22.072478]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:22.073373]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.074091]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.075368]: Fitting sample: DO1004
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:22.077175]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.078749]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:22.084581]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.085791]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:22.087099]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.08794]: 0.024 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.089501]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:22.090224]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:23.047203]: Processing sample `DO1004`.
#> ✔ [2024-09-23 10:00:24.348258]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.349258]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.349986]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:24.350643]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:24.352369]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.353047]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:24.353721]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.354369]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:24.359308]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.360284]: Total 2.3 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Copy Number Alterations (CN48) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.362272]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.363012]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.363738]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.364454]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.366276]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.36776]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:24.368793]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:24.369682]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:24.370596]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.371344]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:24.372082]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.372806]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:24.373592]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:24.375406]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:24.376734]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.37777]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.378764]: Fitting sample: DO1004
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:24.37974]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.3806]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:24.386383]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.3873]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:24.388005]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.388699]: 0.022 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.390089]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:24.390806]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:25.330145]: Processing sample `DO1004`.
#> ✔ [2024-09-23 10:00:26.64639]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:26.647328]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:26.648065]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:26.64869]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:26.650378]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:26.651047]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:26.651652]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:26.65223]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-23 10:00:26.656681]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-23 10:00:26.657545]: Total 2.295 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Write Output ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Raw counts (tally) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ SBS96 tally written to csv: 'pcawg_signature_results/DO1004/SBS96_catalogue.DO1004.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ SBS1536 tally written to csv: 'pcawg_signature_results/DO1004/SBS1536_catalogue.DO1004.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ ID83 tally written to csv: 'pcawg_signature_results/DO1004/ID83_catalogue.DO1004.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ DBS78 tally written to csv: 'pcawg_signature_results/DO1004/DBS78_catalogue.DO1004.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ CN48 tally written to csv: 'pcawg_signature_results/DO1004/CN48_catalogue.DO1004.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ CN80 tally written to csv: 'pcawg_signature_results/DO1004/CN80_catalogue.DO1004.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ CN176 tally written to csv: 'pcawg_signature_results/DO1004/CN176_catalogue.DO1004.hg19.tally.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Fit (Exposures) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ SBS96 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1004/SBS96_fit.DO1004.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ SBS96 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1004/SBS96_fit.DO1004.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ SBS96 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1004/SBS96_fit.DO1004.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ SBS96 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1004/SBS96_fit.DO1004.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ SBS96 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1004/SBS96_fit.DO1004.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ SBS96 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1004/SBS96_fit.DO1004.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ DBS78 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1004/DBS78_fit.DO1004.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ DBS78 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1004/DBS78_fit.DO1004.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ DBS78 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1004/DBS78_fit.DO1004.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ DBS78 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1004/DBS78_fit.DO1004.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ DBS78 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1004/DBS78_fit.DO1004.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ DBS78 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1004/DBS78_fit.DO1004.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ ID83 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1004/ID83_fit.DO1004.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ ID83 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1004/ID83_fit.DO1004.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ ID83 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1004/ID83_fit.DO1004.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ ID83 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1004/ID83_fit.DO1004.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ ID83 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1004/ID83_fit.DO1004.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ ID83 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1004/ID83_fit.DO1004.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ CN48 model fit [expo] has been written to csv: 'pcawg_signature_results/DO1004/CN48_fit.DO1004.hg19.expo.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ CN48 model fit [expo_bootstraps] has been written to csv: 'pcawg_signature_results/DO1004/CN48_fit.DO1004.hg19.expo_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ CN48 model fit [bootstrap_summary] has been written to csv: 'pcawg_signature_results/DO1004/CN48_fit.DO1004.hg19.bootstrap_summary.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ CN48 model fit [error_and_cosine] has been written to csv: 'pcawg_signature_results/DO1004/CN48_fit.DO1004.hg19.error_and_cosine.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ CN48 model fit [error_and_cosine_bootstraps] has been written to csv: 'pcawg_signature_results/DO1004/CN48_fit.DO1004.hg19.error_and_cosine_bootstraps.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ CN48 model fit [p_val] has been written to csv: 'pcawg_signature_results/DO1004/CN48_fit.DO1004.hg19.p_val.csv.gz'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#>  Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100% Progress: ──────────────────────────────────────────────────── 100%
#> ✔ Creating Output Directory at 'pcawg_signature_results/DO1004' [12s]
#> 
#> ℹ Finished successfully
#> ✔ Finished successfully [7ms]
#> 
#> ℹ Cohort analysis completed in 1.08 seconds
#> ✔ Analysis of all 5 samples was successful
```

### Create a reference set

Single sample mutational signature analyses are much more informative
when you have a reference set of samples to compare against.

A reference set allows a single sample analysis to include

1.  Similar sample analysis
2.  UMAPs
3.  

``` r
sig_create_reference_set(
  "pcawg_signature_results", 
  outfolder = "pcawg_reference_set", 
  create_umaps = TRUE, 
  umap_n_neighbours = 2
)
#> 
#> ── Catalogue/Tally Reference Matrix ────────────────────────────────────────────
#> ℹ Reading Catalalogue Tallies
#> ℹ Writing tally refmatrix to pcawg_reference_set/refmatrix.tally.parquet
#> ✔ Tally refmatrix succesfully created
#> 
#> ── UMAPs (from sample catalogues) ──────────────────────────────────────────────
#> ℹ Building UMAPs from sample catalogues
#> 
#> ℹ Building UMAPs from sample catalogues── UMAPs for each signature class ──
#> ℹ Building UMAPs from sample catalogues
#> ℹ Building UMAPs from sample catalogues✔ Building UMAPs from sample catalogues [15ms]
#> 
#> ℹ Bulding UMAP for class CN176
#> Warning: failed creating initial embedding; using random embedding instead
#> Warning: failed creating initial embedding; using random embedding instead
#> ✔ Bulding UMAP for class CN176 [83ms]
#> 
#> ℹ Writing CN176 umap reference to 'pcawg_reference_set/refmatrix.CN176.umap.rds…
#> ✔ Writing CN48 umap reference to 'pcawg_reference_set/refmatrix.CN48.umap.rds.b…
#> 
#> ℹ Bulding UMAP for class CN48
#> Warning: failed creating initial embedding; using random embedding instead
#> Warning: failed creating initial embedding; using random embedding instead
#> ✔ Bulding UMAP for class CN48 [26ms]
#> 
#> ℹ Writing CN48 umap reference to 'pcawg_reference_set/refmatrix.CN48.umap.rds.b…
#> ✔ Writing CN80 umap reference to 'pcawg_reference_set/refmatrix.CN80.umap.rds.b…
#> 
#> ℹ Bulding UMAP for class CN80
#> ✔ Bulding UMAP for class CN80 [27ms]
#> 
#> ℹ Writing CN80 umap reference to 'pcawg_reference_set/refmatrix.CN80.umap.rds.b…
#> ✔ Writing DBS78 umap reference to 'pcawg_reference_set/refmatrix.DBS78.umap.rds…
#> 
#> ℹ Bulding UMAP for class DBS78
#> ✔ Bulding UMAP for class DBS78 [36ms]
#> 
#> ℹ Writing DBS78 umap reference to 'pcawg_reference_set/refmatrix.DBS78.umap.rds…
#> ✔ Writing ID83 umap reference to 'pcawg_reference_set/refmatrix.ID83.umap.rds.b…
#> 
#> ℹ Bulding UMAP for class ID83
#> Warning: failed creating initial embedding; using random embedding instead
#> Warning: failed creating initial embedding; using random embedding instead
#> ✔ Bulding UMAP for class ID83 [26ms]
#> 
#> ℹ Writing ID83 umap reference to 'pcawg_reference_set/refmatrix.ID83.umap.rds.b…
#> ✔ Writing SBS1536 umap reference to 'pcawg_reference_set/refmatrix.SBS1536.umap…
#> 
#> ℹ Bulding UMAP for class SBS1536
#> ✔ Bulding UMAP for class SBS1536 [60ms]
#> 
#> ℹ Writing SBS1536 umap reference to 'pcawg_reference_set/refmatrix.SBS1536.umap…
#> ✔ Writing SBS96 umap reference to 'pcawg_reference_set/refmatrix.SBS96.umap.rds…
#> 
#> ℹ Bulding UMAP for class SBS96
#> ✔ Bulding UMAP for class SBS96 [27ms]
#> 
#> ℹ Writing SBS96 umap reference to 'pcawg_reference_set/refmatrix.SBS96.umap.rds…
#> ── UMAP from all signature classes ──
#> ℹ Writing SBS96 umap reference to 'pcawg_reference_set/refmatrix.SBS96.umap.rds…
#> ℹ Writing SBS96 umap reference to 'pcawg_reference_set/refmatrix.SBS96.umap.rds…✔ Writing SBS96 umap reference to 'pcawg_reference_set/refmatrix.SBS96.umap.rds…
#> 
#> ℹ Building a single UMAP using ALL features accross all signature classes
#> ℹ 5/5 have had features extracted and counted for All signature analyses. Umap will be built using only these 5
#> ℹ Building a single UMAP using ALL features accross all signature classes✔ Building a single UMAP using ALL features accross all signature classes [84ms]
#> 
#> ℹ Writing pan-class umap reference to pcawg_reference_set/refmatrix.panclass.um…
#> ✔ Writing pan-class umap reference to pcawg_reference_set/refmatrix.panclass.um…
#> 
#> 
#> ── Exposure Reference Matrix ───────────────────────────────────────────────────
#> ℹ Reading Signature Exposure Results
#> ✔ Exposure refmatrix successfully created
#> ℹ Reading Signature Exposure Results
#> ✔ Exposure refmatrix successfully created
```
