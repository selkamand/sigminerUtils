devtools::load_all();
noindel = system.file("test.noindels.maf", package = "sigminerUtils")

sig_analyse_mutations(
  maf = noindel,
  ref = "hg19"
)
