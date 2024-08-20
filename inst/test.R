devtools::load_all();
noindel = system.file("test.noindels.maf", package = "sigminerUtils")
nodoublets = system.file("test.nodoublets.maf", package = "sigminerUtils")

sig_analyse_mutations(
  maf = noindel,
  ref = "hg19"
)

sig_analyse_mutations(
  maf = nodoublets,
  ref = "hg19"
)

