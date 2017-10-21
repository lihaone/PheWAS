# script to calculate ePheWAS p values

# input:

# pheno:      an m by n numeric data.frame containing m samples in rows and n phenotypes in columns
# expression: an m by q numeric data.frame of expressions of the gene-of-interest
#             containing m samples in rows and gene expressions in q tissues in columns
# NOTE: pheno and expression data.frame should have the same rownames (samples)
# kmatrix:    an m by m kinship matrix calculated from calculate.kinship function

# output:
# pmatrix.plot: a logpmatrix followed by the Pheno_aligner, for ePheWAS.plot function

require(emma)
ePheWAS.calc <- function(pheno, expression, kmatrix = NULL){
  # check whether pheno and expression have the same rownames
  if(!identical(rownames(pheno), rownames(expression))){
    stop("pheno and expression should have the same rownames")
  }

  # use mixed model to correct for population structure
  if(is.null(kmatrix)){
    stop("kmatrix is needed for emma analysis")
  }else{
    ys <- t(as.matrix(pheno))
    xs <- t(as.matrix(expression))

    # normalize the expression matrix to the range of (0,1)
    range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
    xs <- apply(xs, 1, function(x) range01(x, na.rm = TRUE) )
    xs <- t(xs)

    pmatrix <- emma.ML.LRT(ys=ys, xs=xs, K=kmatrix, ponly=TRUE)
    pmatrix <- data.frame(t(pmatrix), stringsAsFactors = FALSE)
    rownames(pmatrix) <- rownames(ys)
    colnames(pmatrix) <- rownames(xs)
  }
  return(pmatrix)
}
