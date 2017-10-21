# Utility functions for mediation

# prepare for the probability matrix used in mediation

# simple function to convert genotype code into allele probabilities
# input:
# geno.data: a numeric value of the genotype code of a SNP in one sample
single.snp.prob <- function(geno.data){
  # genotype code: -1 - BB, 1 - DD, 0 - BD, 9 - NA
  # return the probability of B and D alleles
  if (!is.numeric(geno.data)){
    stop("geno.data should be a numeric value")
  }
  if (is.na(geno.data)) {
    return(array(data = c(NA, NA)))
  }
  if (geno.data == -1) { # BB
    return(array(data = c(1, 0)))
  }
  if (geno.data == 1){ # DD
    return(array(data = c(0, 1)))
  }
  if (geno.data == 0){ # BD
    return(array(data = c(0.5, 0.5)))
  }
  if (geno.data == 9){ # NA
    return(array(data = c(NA, NA)))
  }
  stop("geno.data should be coded as \"-1, 1, 0, 9\"")
}


# function to convert genotype matrix into genotype probability 3D matrix
# geno:  a p by m numeric genotype matrix containing p markers in rows and m samples in columns,
#     colnames(genotype) should be the sample IDs
#     NOTE: the input genotype matrix should be coded to represent the alleles,
#     i.e. -1 = BB, 1 = DD, 0 = BD, 9 = NA
snp.prob.matrix <- function(geno){

  prob.matrix <- array(NA, c(nrow(geno), ncol(geno), 2),
                dimnames = list(rownames(geno), colnames(geno), c("B", "D")))

  for(row.i in 1:nrow(geno)) {
    for(col.i in 1:ncol(geno)) {
      prob <- single.snp.prob(geno[row.i, col.i])
      prob.matrix[row.i, col.i, 1:2] <- prob
    }
  }
  return(prob.matrix)
}
