devtools::load_all();

path_snvs <- system.file("colo829_testfiles/COLO829v003T.purple.somatic.vcf.gz", package = "sigminerUtils")
path_cnvs <- system.file("colo829_testfiles/COLO829v003T.purple.cnv.somatic.tsv", package = "sigminerUtils")
path_svs <-  system.file("colo829_testfiles/COLO829v003T.purple.sv.vcf.gz", package = "sigminerUtils")

sig_analyse_mutations_single_sample_from_files(
  sample_id = "COLO829v003T", vcf_snv = path_snvs, segment = path_cnvs, vcf_sv = path_svs,
  pass_only = TRUE,
  ref = "hg38",
  output_dir = "colo829_signature_results"
)
