# function to normalize phenotype data

# input:

# pheno.data: numeric phenotype values

pheno.norm <- function(pheno.data){
  # test if there is any zero or negative values in the phenotype data
  if(length(which(pheno.data <= 0)) != 0){
    # if non-positive value exists, use quantile transformation
    message(paste0("Phenotype has non-positive values, quantile transformation will be applied!"))

    quantNorm = function(x){
      qnorm(rank(x, na.last = "keep",ties.method = "average")/(length(x)+1))
    }

    pheno.data.norm <- quantNorm(pheno.data)

  }else{
    # if all values are positive, boxcox transformation could be used
    message(paste0("Phenotype has only positive values, Boxcox transformation will be applied!"))

    # scale all data into data that centers around 1, by dividing the average
    pheno.center <- pheno.data / mean(pheno.data, na.rm = TRUE)

    # run the box-cox transformation
    bc <- boxcox(pheno.center ~ 1, lambda = seq(-2, 2, 0.25), plotit=FALSE)
    trans <- bc$x[which.max(bc$y)]

    if(trans == 0){
      pheno.data.norm <- log(pheno.center)
    }else{
      pheno.data.norm <- (pheno.center^trans - 1)/trans
    }
  }
  return(pheno.data.norm)
}
