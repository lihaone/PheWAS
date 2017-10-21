# script to calculate PheWAS p values

### input:

# pheno: an m by n numeric data.frame containing m samples in rows and n phenotypes in columns
# geno:  a p by m numeric genotype matrix containing p markers in rows and m samples in columns
#     colnames(genotype) should be the sample IDs
#     NOTE: the input genotype matrix should be coded in a minor allele fashion,
#     i.e. minor allel = 0, major allele = 1, heterozygous = 0.5
# kmatrix: an m by m kinship matrix calculated from calculate.kinship function

### output:
# pmatrix containing the PheWAS p values

require(emma)
require(mlmm)
require(plyr)

phewas.calc <- function(pheno, geno, kmatrix, pmatrix=NULL){
  pheno.num <- ncol(pheno)

  for(i in 1:pheno.num){
    pheno.i <- pheno[, i]
    pheno.i.not_na <- !is.na(pheno.i)
    pheno.i <- pheno.i[pheno.i.not_na]
    pheno.length <- length(pheno.i)
    geno.at_index  <- t(geno)[pheno.i.not_na, ]

    if(pheno.length < 15){
      pheno.name <- colnames(pheno)[i]
      message(paste0("phenotype ", pheno.name, " has ", pheno.length, " strains, not enough for analysis!"))
      next
    }else{
      geno_imp.at_index <- as.matrix(geno.at_index)
      kmatrix.i <- kmatrix[pheno.i.not_na, pheno.i.not_na]

      mygwas <- mlmm(Y = pheno.i, X = geno_imp.at_index, K = kmatrix.i, nbchunks = 2, maxsteps = 3)

      pval <- mygwas$pval_step[[1]]$out
      colnames(pval)[2] <- pheno.name
      if (length(pmatrix) == 0){
        pmatrix <- pval
      }else{
        pmatrix <- join(pmatrix, pval, by = "SNP", type = "left", match = "all")
      }
    }
  }
  return(pmatrix)
}
