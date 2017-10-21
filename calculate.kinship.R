# function to calculate kinship matrix from genotype data

# input:

# geno: a p by m numeric genotype matrix containing p markers in rows and m samples in columns
#     colnames(genotype) should be the sample IDs
#     NOTE: the input genotype matrix should be coded in a minor allele fashion,
#     i.e. minor allel = 0, major allele = 1, heterozygous = 0.5

require(emma)
kinship.emma <- function(geno) {
  K_mat <- emma.kinship(geno)            # Calculate kinship using the emma package
  colnames(K_mat) <- colnames(geno)      # emma does not set the rownames and column names, so we should do it
  rownames(K_mat) <- colnames(geno)
  return(K_mat)
}
