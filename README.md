
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
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:48.084332]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:50.951023]: We would assume you marked all variants' position in + strand.
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
#> ✔ [2024-09-20 10:42:51.022241]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.024789]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.026161]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.028527]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.029251]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.031549]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.03318]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.034245]: SBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.03525]: Extracting 5' and 3' adjacent bases.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.467732]: Extracting +/- 20bp around mutated bases for background C>T estimation.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.669191]: Estimating APOBEC enrichment scores.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.67009]: Performing one-way Fisher's test for APOBEC enrichment.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.672491]: APOBEC related mutations are enriched in 0% of samples (APOBEC enrichment score > 2; 0 of 1 samples)
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.673856]: Creating SBS sample-by-component matrices.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.675727]: SBS-6 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.67746]: SBS-96 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.682764]: SBS-1536 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.683374]: Return SBS-96 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.684437]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.685194]: 3.601 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.685825]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.689734]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.690425]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.692278]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.693687]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.69599]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.696679]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.698787]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.700403]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.701178]: DBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.701882]: Searching DBS records...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.711654]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.809051]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.813405]: DBS-78 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.817806]: DBS-1248 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.818501]: Return SBS-78 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.819439]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.820038]: 0.134 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.820648]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:51.823761]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.824497]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.82631]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.827548]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.829685]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.830354]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:51.832253]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.098852]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.0999]: INDEL matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.300723]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.301774]: INDEL length extracted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.302863]: Adjacent copies counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.304394]: Microhomology size calculated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.308484]: INDEL records classified into different components (types).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.310489]: ID-28 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.312284]: ID-83 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.312893]: Return ID-83 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.313561]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.314173]: 0.494 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Tally (CopyNumber) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.316072]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.316675]: Genome build  : hg19.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.317254]: Genome measure: wg.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.317838]: When add_loh is TRUE, use_all is forced to TRUE.
#> Please drop columns you don't want to keep before reading.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'sigminer version 2.3.2
#> - Star me at https://github.com/ShixiangWang/sigminer
#> - Run hello() to see usage and citation.
#> ✔ [2024-09-20 10:42:52.323336]: Chromosome size database for build obtained.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.324033]: Reading input.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.324631]: A data frame as input detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.325248]: Column names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.325883]: Column order set.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.326612]: Rows with NA copy number removed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.327932]: Chromosomes unified.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.329235]: Value 2 (normal copy) filled to uncalled chromosomes.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.330678]: Data imported.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.331294]: Segments info:
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.331879]:     Keep - 253
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.332483]:     Drop - 0
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.333278]: Segments sorted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.333864]: Adding LOH labels...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.334653]: Skipped joining adjacent segments with same copy number value.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.335361]: Segmental table cleaned.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.335987]: Annotating.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.341569]: Annotation done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.342209]: Summarizing per sample.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.348333]: Summarized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.348999]: Generating CopyNumber object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.349871]: Generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.35047]: Validating object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.351086]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.351711]: 0.036 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.35232]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.35363]: When you use method 'S', please make sure you have set 'join_adj_seg' to FALSE and 'add_loh' to TRUE in 'read_copynumber() in the previous step!
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.360335]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.361106]: 0.009 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.36173]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.362645]: Step: getting copy number features.
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
#> ℹ [2024-09-20 10:42:52.396608]: Getting breakpoint count per 10 Mb...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.424923]: Getting breakpoint count per chromosome arm...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.432738]: Getting copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.4346]: Getting change-point copy number change...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.441123]: Getting length of chains of oscillating copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.448019]: Getting (log10 based) segment size...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.449789]: Getting the minimal number of chromosome with 50% CNV...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.45445]: Getting burden of chromosome...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.463085]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.46921]: Step: generating copy number components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.469808]: `feature_setting` checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.471068]: Step: counting components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.558411]: Counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.55927]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.560118]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.560713]: 0.199 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.561286]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.561858]: Generated 'Index' column to track the copy number segment location.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.562753]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.564864]: NOTE: this method derives features for each segment. Be patient...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.571608]: Getting absolute copy number value of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.578309]: Getting segment size of eash segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.579863]: Getting context change shape based on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.58691]: Getting change extent on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.593837]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.594501]: Step: generating copy number components based on combination.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.595551]: Classified and combined.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.596114]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.60906]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.609716]: 0.048 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Fitting ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Single Base Substitutions (SBS96) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.612984]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.61354]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.614103]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.614676]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.61622]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.616764]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.617321]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.61803]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.618906]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.619489]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.620068]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.620637]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.621184]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.621733]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.626744]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.627336]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.628097]: Fitting sample: DO1000
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.630382]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.631052]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.634778]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.635369]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:52.636042]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.636646]: 0.02 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.637847]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.638408]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:52.648268]: Processing sample `DO1000`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.362932]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.363695]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.364369]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.364943]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.36881]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.369467]: 0.005 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.370048]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.370604]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.375762]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.376694]: Total 1.764 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── INDELS (ID83) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.378054]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.378628]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.379198]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.379776]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.381104]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.381702]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.382278]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.3829]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.383566]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.384143]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.384705]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.385265]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.385826]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.38639]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.386957]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.387515]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.388166]: Fitting sample: DO1000
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.388826]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.389398]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.392744]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.393344]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:54.393964]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.394546]: 0.013 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.395677]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:54.396255]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:55.719223]: Processing sample `DO1000`.
#> ✔ [2024-09-20 10:42:57.017177]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.017995]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.018617]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:57.019166]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:57.020574]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.021146]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:57.021712]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.022242]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:57.026806]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.027703]: Total 2.65 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Doublet Mutations (DBS78) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.030211]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.033192]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.035144]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.035978]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.037516]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.03815]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:57.038802]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:57.039493]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:57.040249]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.041741]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:57.042536]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.043192]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:57.044274]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:57.045001]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:57.045654]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.046762]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.047615]: Fitting sample: DO1000
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:57.049476]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.051083]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:57.055085]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.055787]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:57.056447]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.05705]: 0.02 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.058261]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:57.058846]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:58.004806]: Processing sample `DO1000`.
#> ✔ [2024-09-20 10:42:59.308671]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.309479]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.310102]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:59.310648]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:59.311956]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.312522]: 0.002 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:59.313084]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.313618]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:59.317944]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.318773]: Total 2.289 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'── Copy Number Alterations (CN48) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.320481]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.321122]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.321743]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.322353]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.323948]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.324611]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:59.325292]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:59.326042]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:59.326822]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.327486]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:59.328118]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.328751]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:59.329365]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:59.329983]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:59.330613]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.331241]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.332005]: Fitting sample: DO1000
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:59.332933]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.334029]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:59.338998]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.34063]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:42:59.341429]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.342135]: 0.018 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.344881]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:42:59.345576]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:43:00.278586]: Processing sample `DO1000`.
#> ✔ [2024-09-20 10:43:01.592857]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:43:01.593739]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:43:01.594446]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:43:01.595097]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:43:01.596809]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:43:01.597469]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:43:01.598138]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:43:01.598762]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'✔ [2024-09-20 10:43:01.603623]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1000'ℹ [2024-09-20 10:43:01.604564]: Total 2.284 secs elapsed.
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
#>  Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100%
#> ✔ Creating Output Directory at 'pcawg_signature_results/DO1000' [15.5s]
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
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:02.805077]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:02.80931]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:02.810108]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:02.812548]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:02.815053]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:02.820409]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:02.82122]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:02.825785]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:02.827657]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:02.828261]: SBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:02.829157]: Extracting 5' and 3' adjacent bases.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.436905]: Extracting +/- 20bp around mutated bases for background C>T estimation.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.607496]: Estimating APOBEC enrichment scores.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.608403]: Performing one-way Fisher's test for APOBEC enrichment.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.610459]: APOBEC related mutations are enriched in 0% of samples (APOBEC enrichment score > 2; 0 of 1 samples)
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.611733]: Creating SBS sample-by-component matrices.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.613507]: SBS-6 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.615507]: SBS-96 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.627756]: SBS-1536 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.628595]: Return SBS-96 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.62965]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.630235]: 0.825 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.630807]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.634204]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.634966]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.637212]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.639656]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.646208]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.646953]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.651878]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.653518]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.654083]: DBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.654836]: Searching DBS records...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.676042]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.831621]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.836001]: DBS-78 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.840493]: DBS-1248 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.841138]: Return SBS-78 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.842033]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.842622]: 0.212 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.843199]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.846457]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.847184]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.849467]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.851987]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.857853]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.858678]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.863654]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:03.865378]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:03.865945]: INDEL matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.152263]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.153285]: INDEL length extracted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.154873]: Adjacent copies counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.166508]: Microhomology size calculated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.172892]: INDEL records classified into different components (types).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.174593]: ID-28 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.176264]: ID-83 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.176839]: Return ID-83 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.177489]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.178098]: 0.335 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Tally (CopyNumber) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.179586]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.180149]: Genome build  : hg19.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.180697]: Genome measure: wg.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.181287]: When add_loh is TRUE, use_all is forced to TRUE.
#> Please drop columns you don't want to keep before reading.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.182951]: Chromosome size database for build obtained.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.183537]: Reading input.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.184146]: A data frame as input detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.184742]: Column names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.185336]: Column order set.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.186017]: Rows with NA copy number removed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.187673]: Chromosomes unified.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.18895]: Value 2 (normal copy) filled to uncalled chromosomes.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.190287]: Data imported.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.190877]: Segments info:
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.191434]:     Keep - 637
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.19202]:     Drop - 0
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.192768]: Segments sorted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.193315]: Adding LOH labels...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.19405]: Skipped joining adjacent segments with same copy number value.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.194678]: Segmental table cleaned.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.195269]: Annotating.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.200984]: Annotation done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.201626]: Summarizing per sample.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.206071]: Summarized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.206691]: Generating CopyNumber object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.207381]: Generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.20798]: Validating object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.208581]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.209188]: 0.03 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.209756]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.210627]: When you use method 'S', please make sure you have set 'join_adj_seg' to FALSE and 'add_loh' to TRUE in 'read_copynumber() in the previous step!
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.215606]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.216238]: 0.006 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.216816]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.217684]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.230143]: Getting breakpoint count per 10 Mb...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.244886]: Getting breakpoint count per chromosome arm...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.25228]: Getting copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.253906]: Getting change-point copy number change...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.25998]: Getting length of chains of oscillating copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.267118]: Getting (log10 based) segment size...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.268863]: Getting the minimal number of chromosome with 50% CNV...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.272766]: Getting burden of chromosome...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.276286]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.276904]: Step: generating copy number components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.27747]: `feature_setting` checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.278783]: Step: counting components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.374373]: Counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.375261]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.376139]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.376759]: 0.16 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.377335]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.377935]: Generated 'Index' column to track the copy number segment location.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.37889]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.381015]: NOTE: this method derives features for each segment. Be patient...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.388952]: Getting absolute copy number value of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.395888]: Getting segment size of eash segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.397481]: Getting context change shape based on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.404863]: Getting change extent on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.412059]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.412679]: Step: generating copy number components based on combination.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.41362]: Classified and combined.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.414178]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.426444]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.427165]: 0.05 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Fitting ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Single Base Substitutions (SBS96) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.430289]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.430857]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.431419]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.432019]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.433392]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.433982]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.43458]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.435298]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.436128]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.436716]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.437308]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.437885]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.438461]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.439037]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.439614]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.440166]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.440793]: Fitting sample: DO1001
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.442574]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.443112]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.446419]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.447034]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:04.44767]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.448251]: 0.015 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.449362]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.449953]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:04.46032]: Processing sample `DO1001`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.163534]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.164494]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.165139]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.165689]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.169124]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.169739]: 0.005 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.170347]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.170893]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.175837]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.176727]: Total 1.746 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── INDELS (ID83) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.178049]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.178639]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.179223]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.179816]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.181119]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.181677]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.182273]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.182915]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.183608]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.184177]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.184714]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.185252]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.185799]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.186348]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.186915]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.187492]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.18814]: Fitting sample: DO1001
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.188819]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.189363]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.192531]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.193093]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:06.193701]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.194292]: 0.013 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.195423]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:06.195975]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:07.157845]: Processing sample `DO1001`.
#> ✔ [2024-09-20 10:43:08.487868]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.488701]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.489335]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:08.489903]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:08.491313]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.491934]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:08.492541]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.493137]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:08.49749]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.498322]: Total 2.32 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Doublet Mutations (DBS78) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.500049]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.500695]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.501822]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.504584]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.507439]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.508175]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:08.508842]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:08.509531]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:08.510271]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.510899]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:08.511569]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.51227]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:08.512936]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:08.51382]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:08.515267]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.516008]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.516831]: Fitting sample: DO1001
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:08.517598]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.51825]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:08.524778]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.525744]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:08.526461]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.527065]: 0.02 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.52835]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:08.528962]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:09.482505]: Processing sample `DO1001`.
#> ✔ [2024-09-20 10:43:10.875809]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.876701]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.87741]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:10.878056]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:10.879576]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.880243]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:10.880878]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.881513]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:10.886273]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.887155]: Total 2.387 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'── Copy Number Alterations (CN48) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.888995]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.889686]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.890358]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.89236]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.896014]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.896849]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:10.897618]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:10.898442]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:10.899321]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.900078]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:10.900826]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.901554]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:10.902916]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:10.904065]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:10.90486]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.905748]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.90757]: Fitting sample: DO1001
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:10.908537]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.909689]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:10.914922]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.915744]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:10.916403]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.917044]: 0.021 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.918346]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:10.918945]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:11.884995]: Processing sample `DO1001`.
#> ✔ [2024-09-20 10:43:13.225005]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:13.226024]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:13.22666]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:13.227214]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:13.228725]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:13.229318]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:13.229916]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:13.230493]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'✔ [2024-09-20 10:43:13.234856]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1001'ℹ [2024-09-20 10:43:13.235683]: Total 2.347 secs elapsed.
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
#>  Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100%
#> ✔ Creating Output Directory at 'pcawg_signature_results/DO1001' [11.9s]
#> 
#> ℹ Finished successfully
#> ✔ Finished successfully [6ms]
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
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:14.628444]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:14.632241]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:14.633041]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:14.635549]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:14.638025]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:14.643311]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:14.644125]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:14.648727]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:14.65054]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:14.651145]: SBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:14.652034]: Extracting 5' and 3' adjacent bases.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.24579]: Extracting +/- 20bp around mutated bases for background C>T estimation.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.42011]: Estimating APOBEC enrichment scores.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.421019]: Performing one-way Fisher's test for APOBEC enrichment.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.423127]: APOBEC related mutations are enriched in 0% of samples (APOBEC enrichment score > 2; 0 of 1 samples)
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.424382]: Creating SBS sample-by-component matrices.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.426114]: SBS-6 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.427837]: SBS-96 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.433281]: SBS-1536 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.433925]: Return SBS-96 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.435004]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.43561]: 0.807 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.436208]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.439372]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.440116]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.442479]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.445113]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.450831]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.451652]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.456277]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.457934]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.458499]: DBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.459321]: Searching DBS records...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.480113]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.705919]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.71028]: DBS-78 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.720573]: DBS-1248 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.721423]: Return SBS-78 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.722334]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.722933]: 0.287 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.723582]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.726822]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.727586]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.729956]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.732598]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.739272]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.740091]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.744759]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:15.746561]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:15.747125]: INDEL matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.056411]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.057592]: INDEL length extracted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.059874]: Adjacent copies counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.086514]: Microhomology size calculated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.096363]: INDEL records classified into different components (types).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.098298]: ID-28 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.100107]: ID-83 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.100734]: Return ID-83 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.101403]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.102036]: 0.378 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Tally (CopyNumber) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.103597]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.10421]: Genome build  : hg19.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.104806]: Genome measure: wg.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.105395]: When add_loh is TRUE, use_all is forced to TRUE.
#> Please drop columns you don't want to keep before reading.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.107039]: Chromosome size database for build obtained.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.107655]: Reading input.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.108267]: A data frame as input detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.108904]: Column names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.109534]: Column order set.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.110192]: Rows with NA copy number removed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.111755]: Chromosomes unified.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.113028]: Value 2 (normal copy) filled to uncalled chromosomes.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.114464]: Data imported.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.115075]: Segments info:
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.115671]:     Keep - 573
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.116268]:     Drop - 0
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.116988]: Segments sorted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.117578]: Adding LOH labels...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.118362]: Skipped joining adjacent segments with same copy number value.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.119]: Segmental table cleaned.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.119588]: Annotating.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.125276]: Annotation done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.125917]: Summarizing per sample.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.13054]: Summarized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.131194]: Generating CopyNumber object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.131921]: Generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.132549]: Validating object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.133157]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.133767]: 0.03 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.134378]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.135327]: When you use method 'S', please make sure you have set 'join_adj_seg' to FALSE and 'add_loh' to TRUE in 'read_copynumber() in the previous step!
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.140576]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.141238]: 0.007 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.141856]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.142754]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.155928]: Getting breakpoint count per 10 Mb...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.170887]: Getting breakpoint count per chromosome arm...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.178309]: Getting copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.180017]: Getting change-point copy number change...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.186272]: Getting length of chains of oscillating copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.193168]: Getting (log10 based) segment size...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.194863]: Getting the minimal number of chromosome with 50% CNV...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.198873]: Getting burden of chromosome...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.202656]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.20335]: Step: generating copy number components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.203942]: `feature_setting` checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.205233]: Step: counting components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.293584]: Counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.294578]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.295418]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.296035]: 0.154 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.296659]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.297301]: Generated 'Index' column to track the copy number segment location.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.298285]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.300365]: NOTE: this method derives features for each segment. Be patient...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.308526]: Getting absolute copy number value of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.31572]: Getting segment size of eash segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.317342]: Getting context change shape based on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.324855]: Getting change extent on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.332315]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.332983]: Step: generating copy number components based on combination.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.333969]: Classified and combined.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.33458]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.347374]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.3482]: 0.052 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Fitting ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Single Base Substitutions (SBS96) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.351459]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.352083]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.352672]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.353245]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.354668]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.355281]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.3559]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.356693]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.357566]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.35816]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.35878]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.35938]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.359988]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.360598]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.361219]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.361785]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.362463]: Fitting sample: DO1002
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.364354]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.364916]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.368568]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.369209]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:16.369877]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.370519]: 0.016 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.371766]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.372357]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:16.383222]: Processing sample `DO1002`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.12822]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.129134]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.12986]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.130477]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.134145]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.134861]: 0.005 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.135476]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.136063]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.1412]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.142151]: Total 1.791 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── INDELS (ID83) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.143586]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.144186]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.144789]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.145379]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.146738]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.147328]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.14793]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.148595]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.149241]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.149815]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.150402]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.150988]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.151572]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.152153]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.152736]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.153318]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.153996]: Fitting sample: DO1002
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.154697]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.155266]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.158602]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.159212]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:18.159858]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.160468]: 0.014 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.161662]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:18.162281]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:19.166745]: Processing sample `DO1002`.
#> ✔ [2024-09-20 10:43:20.553256]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.554229]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.554877]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:20.555445]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:20.556891]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.557508]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:20.558128]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.558729]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:20.563086]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.563941]: Total 2.42 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Doublet Mutations (DBS78) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.567376]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.57037]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.571479]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.57245]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.574033]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.574691]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:20.575362]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:20.576098]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:20.576879]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.577904]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:20.578592]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.579248]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:20.580551]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:20.581424]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:20.582147]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.582834]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.583659]: Fitting sample: DO1002
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:20.585108]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.587002]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:20.591425]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.592256]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:20.592947]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.593572]: 0.02 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.594818]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:20.595432]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:21.556624]: Processing sample `DO1002`.
#> ✔ [2024-09-20 10:43:22.897105]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.897994]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.898744]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:22.89938]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:22.900856]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.901513]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:22.902153]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.902769]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:22.907558]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.908472]: Total 2.341 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'── Copy Number Alterations (CN48) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.910456]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.912319]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.913134]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.914001]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.915969]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.916796]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:22.917582]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:22.918378]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:22.919197]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.919907]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:22.920635]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.9214]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:22.922194]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:22.922917]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:22.923859]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.925153]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.926066]: Fitting sample: DO1002
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:22.926944]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.92875]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:22.936322]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.937309]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:22.938014]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.938692]: 0.023 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.940094]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:22.94075]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:23.90787]: Processing sample `DO1002`.
#> ✔ [2024-09-20 10:43:25.267662]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:25.2686]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:25.269251]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:25.269826]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:25.271314]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:25.27193]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:25.272556]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:25.273159]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'✔ [2024-09-20 10:43:25.277621]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1002'ℹ [2024-09-20 10:43:25.278462]: Total 2.368 secs elapsed.
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
#>  Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100%
#> ✔ Creating Output Directory at 'pcawg_signature_results/DO1002' [11.8s]
#> 
#> ℹ Finished successfully
#> ✔ Finished successfully [6ms]
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
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:26.373191]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:26.377022]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:26.377856]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:26.380206]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:26.382432]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:26.387021]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:26.387836]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:26.391993]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:26.393834]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:26.39448]: SBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:26.395374]: Extracting 5' and 3' adjacent bases.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:26.943705]: Extracting +/- 20bp around mutated bases for background C>T estimation.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.108069]: Estimating APOBEC enrichment scores.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.109056]: Performing one-way Fisher's test for APOBEC enrichment.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.111134]: APOBEC related mutations are enriched in 0% of samples (APOBEC enrichment score > 2; 0 of 1 samples)
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.11243]: Creating SBS sample-by-component matrices.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.114163]: SBS-6 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.116047]: SBS-96 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.121706]: SBS-1536 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.122401]: Return SBS-96 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.123494]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.124165]: 0.751 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.124849]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.128253]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.129067]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.131564]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.133834]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.138604]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.139426]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.143562]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.145294]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.145934]: DBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.146754]: Searching DBS records...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.168877]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.368568]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.373106]: DBS-78 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.377748]: DBS-1248 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.378512]: Return SBS-78 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.379513]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.380155]: 0.255 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.380791]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.384007]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.384771]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.387035]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.389199]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.394093]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.394887]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.39902]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.4009]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.401485]: INDEL matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.671929]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.673014]: INDEL length extracted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.674226]: Adjacent copies counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.682745]: Microhomology size calculated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.68841]: INDEL records classified into different components (types).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.690222]: ID-28 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.691937]: ID-83 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.692573]: Return ID-83 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.693259]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.693904]: 0.313 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Tally (CopyNumber) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.695512]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.696137]: Genome build  : hg19.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.696757]: Genome measure: wg.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.69736]: When add_loh is TRUE, use_all is forced to TRUE.
#> Please drop columns you don't want to keep before reading.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.698997]: Chromosome size database for build obtained.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.699584]: Reading input.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.700199]: A data frame as input detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.700863]: Column names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.70151]: Column order set.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.702196]: Rows with NA copy number removed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.704287]: Chromosomes unified.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.705562]: Value 2 (normal copy) filled to uncalled chromosomes.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.70699]: Data imported.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.707662]: Segments info:
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.708238]:     Keep - 1007
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.708845]:     Drop - 0
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.709626]: Segments sorted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.710242]: Adding LOH labels...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.711057]: Skipped joining adjacent segments with same copy number value.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.711728]: Segmental table cleaned.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.71233]: Annotating.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.718745]: Annotation done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.719457]: Summarizing per sample.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.724061]: Summarized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.724715]: Generating CopyNumber object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.725414]: Generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.726029]: Validating object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.726637]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.727255]: 0.032 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.727877]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.728853]: When you use method 'S', please make sure you have set 'join_adj_seg' to FALSE and 'add_loh' to TRUE in 'read_copynumber() in the previous step!
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.734012]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.734677]: 0.007 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.735309]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.736248]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.748694]: Getting breakpoint count per 10 Mb...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.763082]: Getting breakpoint count per chromosome arm...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.770256]: Getting copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.771863]: Getting change-point copy number change...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.778032]: Getting length of chains of oscillating copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.785276]: Getting (log10 based) segment size...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.786911]: Getting the minimal number of chromosome with 50% CNV...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.790784]: Getting burden of chromosome...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.794631]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.795277]: Step: generating copy number components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.795867]: `feature_setting` checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.797253]: Step: counting components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.888415]: Counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.889383]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.890245]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.890912]: 0.156 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.891573]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.892239]: Generated 'Index' column to track the copy number segment location.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.893232]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.895387]: NOTE: this method derives features for each segment. Be patient...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.903367]: Getting absolute copy number value of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.91046]: Getting segment size of eash segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.912131]: Getting context change shape based on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.919772]: Getting change extent on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.927485]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.928147]: Step: generating copy number components based on combination.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.929227]: Classified and combined.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.929876]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.94288]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.943681]: 0.052 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Fitting ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Single Base Substitutions (SBS96) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.947091]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.947741]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.948381]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.948959]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.950341]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.950983]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.951611]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.952408]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.953344]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.953979]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.954605]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.955232]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.955867]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.956488]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.957124]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.957745]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.958458]: Fitting sample: DO1003
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.960336]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.960914]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.964501]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.965217]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:27.965917]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.972807]: 0.022 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.974806]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.975425]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:27.986117]: Processing sample `DO1003`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.723543]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.724372]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.725016]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.725602]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.729254]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.729917]: 0.005 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.730536]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.731132]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.736245]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.737269]: Total 1.79 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── INDELS (ID83) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.738764]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.739402]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.740016]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.740625]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.742004]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.742604]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.743218]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.743885]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.744596]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.745198]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.745801]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.746402]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.747009]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.747631]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.748233]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.748834]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.749513]: Fitting sample: DO1003
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.750213]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.750815]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.754041]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.754651]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:29.755307]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.75593]: 0.014 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.757084]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:29.757723]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:30.718239]: Processing sample `DO1003`.
#> ✔ [2024-09-20 10:43:32.064666]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.065541]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.066196]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:32.066782]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:32.068303]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.068952]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:32.06959]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.070204]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:32.074657]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.075534]: Total 2.337 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Doublet Mutations (DBS78) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.077401]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.078094]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.078758]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.079423]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.08104]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.08176]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:32.08248]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:32.083268]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:32.084086]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.084783]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:32.08546]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.086118]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:32.086794]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:32.087457]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:32.088132]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.088819]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.089732]: Fitting sample: DO1003
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:32.091396]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.092222]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:32.098214]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.099377]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:32.100183]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.100898]: 0.02 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.102229]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:32.102836]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:33.078726]: Processing sample `DO1003`.
#> ✔ [2024-09-20 10:43:34.448684]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.449832]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.450572]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:34.451237]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:34.452766]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.453467]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:34.454134]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.454784]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:34.459726]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.460728]: Total 2.383 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'── Copy Number Alterations (CN48) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.462786]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.463557]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.467176]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.469365]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.4712]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.471932]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:34.472709]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:34.473582]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:34.474479]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.476337]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:34.477663]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.478619]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:34.48012]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:34.481006]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:34.481789]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.482543]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.484668]: Fitting sample: DO1003
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:34.485852]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.486754]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:34.490998]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.491781]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:34.492508]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.493195]: 0.022 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.49455]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:34.495219]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:35.476222]: Processing sample `DO1003`.
#> ✔ [2024-09-20 10:43:36.848182]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:36.849053]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:36.849725]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:36.850313]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:36.851825]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:36.852436]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:36.853027]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:36.85359]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'✔ [2024-09-20 10:43:36.858032]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1003'ℹ [2024-09-20 10:43:36.858902]: Total 2.396 secs elapsed.
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
#>  Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100%
#> ✔ Creating Output Directory at 'pcawg_signature_results/DO1003' [11.6s]
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
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:38.14051]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:38.144351]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:38.145147]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:38.148592]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:38.152925]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:38.162835]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:38.163844]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:38.17295]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:38.174876]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:38.175576]: SBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:38.176801]: Extracting 5' and 3' adjacent bases.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:38.895054]: Extracting +/- 20bp around mutated bases for background C>T estimation.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.120877]: Estimating APOBEC enrichment scores.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.121932]: Performing one-way Fisher's test for APOBEC enrichment.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.124858]: APOBEC related mutations are enriched in 0% of samples (APOBEC enrichment score > 2; 0 of 1 samples)
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.126286]: Creating SBS sample-by-component matrices.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.128162]: SBS-6 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.130155]: SBS-96 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.136084]: SBS-1536 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.136956]: Return SBS-96 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.138064]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.138719]: 0.998 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.139389]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.142716]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.143496]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.153469]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.158332]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.168603]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.169737]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.178822]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.18082]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.181449]: DBS matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.182462]: Searching DBS records...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.226011]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.504915]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.510245]: DBS-78 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.515724]: DBS-1248 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.516645]: Return SBS-78 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.517821]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.518568]: 0.379 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.519345]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.52318]: We would assume you marked all variants' position in + strand.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.52415]: Reference genome loaded.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.527779]: Variants from MAF object queried.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.53251]: Chromosome names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.545716]: Sex chromosomes properly handled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.547034]: Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.557202]: Variant start and end position checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.559414]: Variant data for matrix generation preprocessed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.560117]: INDEL matrix generation - start.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.894637]: Reference sequences queried from genome.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.895817]: INDEL length extracted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.898672]: Adjacent copies counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.927424]: Microhomology size calculated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.946033]: INDEL records classified into different components (types).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.948112]: ID-28 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.949994]: ID-83 matrix created.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.950655]: Return ID-83 as major matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.951357]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.952013]: 0.433 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Tally (CopyNumber) ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.953678]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.95431]: Genome build  : hg19.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.954927]: Genome measure: wg.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.955542]: When add_loh is TRUE, use_all is forced to TRUE.
#> Please drop columns you don't want to keep before reading.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.957216]: Chromosome size database for build obtained.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.957864]: Reading input.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.958501]: A data frame as input detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.959167]: Column names checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.959827]: Column order set.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.960513]: Rows with NA copy number removed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.96197]: Chromosomes unified.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.963264]: Value 2 (normal copy) filled to uncalled chromosomes.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.964708]: Data imported.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.9654]: Segments info:
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.966027]:     Keep - 442
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.966644]:     Drop - 0
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.967396]: Segments sorted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.96801]: Adding LOH labels...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.968827]: Skipped joining adjacent segments with same copy number value.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.969511]: Segmental table cleaned.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.970139]: Annotating.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.975749]: Annotation done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.976423]: Summarizing per sample.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.981059]: Summarized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.981756]: Generating CopyNumber object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.982514]: Generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.98314]: Validating object.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.983764]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.9844]: 0.031 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.985041]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.985993]: When you use method 'S', please make sure you have set 'join_adj_seg' to FALSE and 'add_loh' to TRUE in 'read_copynumber() in the previous step!
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:39.991234]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.991924]: 0.007 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.992577]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:39.993518]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.008194]: Getting breakpoint count per 10 Mb...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.023011]: Getting breakpoint count per chromosome arm...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.030447]: Getting copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.032087]: Getting change-point copy number change...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.038173]: Getting length of chains of oscillating copy number...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.044812]: Getting (log10 based) segment size...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.046469]: Getting the minimal number of chromosome with 50% CNV...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.05038]: Getting burden of chromosome...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.054479]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.055144]: Step: generating copy number components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.055744]: `feature_setting` checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.057152]: Step: counting components.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.148789]: Counted.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.149766]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.150663]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.151387]: 0.159 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.152048]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.152726]: Generated 'Index' column to track the copy number segment location.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.154161]: Step: getting copy number features.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.156337]: NOTE: this method derives features for each segment. Be patient...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.167077]: Getting absolute copy number value of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.174124]: Getting segment size of eash segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.175752]: Getting context change shape based on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.183426]: Getting change extent on left and right sides of each segment...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.191133]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.191809]: Step: generating copy number components based on combination.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.192863]: Classified and combined.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.193528]: Step: generating components by sample matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.20638]: Matrix generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.207156]: 0.055 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Fitting ──
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Single Base Substitutions (SBS96) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.210614]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.211267]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.211903]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.21254]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.21396]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.214593]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.215243]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.216068]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.217011]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.217654]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.218301]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.21893]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.219562]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.220198]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.220832]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.221473]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.222183]: Fitting sample: DO1004
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.224121]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.224715]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.228255]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.228919]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:40.229617]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.230279]: 0.016 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.231506]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.232168]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:40.244931]: Processing sample `DO1004`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.09751]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.098414]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.099106]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.099726]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.103378]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.104059]: 0.005 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.104697]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.105308]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.110649]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.111626]: Total 1.901 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── INDELS (ID83) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.113128]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.113753]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.114372]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.114988]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.1164]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.117014]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.117646]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.118321]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.119038]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.119656]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.120271]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.120884]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.121496]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.122108]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.122717]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.123342]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.124026]: Fitting sample: DO1004
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.124744]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.125365]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.128757]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.129399]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:42.130062]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.13069]: 0.014 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.131856]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:42.132492]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:43.100449]: Processing sample `DO1004`.
#> ✔ [2024-09-20 10:43:44.457561]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.458477]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.459146]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:44.459742]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:44.461318]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.461947]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:44.462571]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.463184]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:44.467714]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.468595]: Total 2.355 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Doublet Mutations (DBS78) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.470474]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.471215]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.471893]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.472571]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.474217]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.474979]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:44.475746]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:44.476553]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:44.477402]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.478092]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:44.479208]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.480056]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:44.4808]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:44.48202]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:44.482722]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.483418]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.484832]: Fitting sample: DO1004
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:44.485913]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.487111]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:44.493698]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.494677]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:44.495416]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.496074]: 0.022 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.497372]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:44.498009]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:45.465041]: Processing sample `DO1004`.
#> ✔ [2024-09-20 10:43:46.810233]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.811185]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.811909]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:46.812583]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:46.814125]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.814818]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:46.815494]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.816141]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:46.820877]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.82181]: Total 2.351 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'── Copy Number Alterations (CN48) 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.825793]: Batch Bootstrap Signature Exposure Analysis Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.828357]: Samples to be filtered out: 
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.829509]: Finding optimal exposures (&errors) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.830372]: Calling method `QP`.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.832034]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.832777]: Signature index not detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:46.833518]: Signature matrix/data.frame detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:46.834602]: Database and index checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:46.835485]: Signature normalized.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.836233]: Checking row number for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:46.837089]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.838565]: Checking rownames for catalog matrix and signature matrix.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:46.839396]: Checked.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:46.84017]: Method 'QP' detected.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:46.840946]: Corresponding function generated.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.841694]: Calling function.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.842964]: Fitting sample: DO1004
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:46.845814]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.847088]: Generating output signature exposures.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:46.851357]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.852199]: Calculating errors (Frobenius Norm).
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:46.85291]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.853564]: 0.022 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.854893]: Getting bootstrap exposures (&errors/similarity) for different methods.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:46.855568]: This step is time consuming, please be patient.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:47.833315]: Processing sample `DO1004`.
#> ✔ [2024-09-20 10:43:49.19382]: Gotten.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:49.194721]: Reporting p values...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:49.195401]: Started.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:49.195995]: Batch mode enabled.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:49.197632]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:49.198303]: 0.003 secs elapsed.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:49.198957]: Done.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:49.199577]: Cleaning results...
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'✔ [2024-09-20 10:43:49.204082]: Outputing.
#> ℹ Creating Output Directory at 'pcawg_signature_results/DO1004'ℹ [2024-09-20 10:43:49.204979]: Total 2.379 secs elapsed.
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
#>  Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100% Progress: ─────────────────────────────────────────────────────────────────── 100%
#> ✔ Creating Output Directory at 'pcawg_signature_results/DO1004' [12.4s]
#> 
#> ℹ Finished successfully
#> ✔ Finished successfully [6ms]
#> 
#> ℹ Cohort analysis completed in 1.05 seconds
#> ✔ Analysis of all 5 samples was successful
```
