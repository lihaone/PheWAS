# function to calculate reverse-mediation results

require(intermediate)
require(plyr)

# target.gene: the gene symbol of the target gene
# qtl.snp: the SNP ID of the cis-QTL of target gene
# geno:  a p by m numeric genotype matrix containing p markers in rows and m samples in columns
#     colnames(genotype) should be the sample IDs
#     NOTE: the input genotype matrix should be coded in a minor allele fashion,
#     i.e. minor allel = 0, major allele = 1, heterozygous = 0.5
# expression: an m by q numeric data.frame of all gene expressions,
#             containing m samples in rows and q genes in columns
# annotation: an q by 3 data.frame containing the genetic positions of all q genes,
#             column 1 (SYMBOL) = gene symbol, column 2 (CHR) = chromosome, column 3 (POS) = position in Mb
# NOTE: colnames(expression) and annotation[,1] should be identical
# covar: a numeric matrix with additive covariates

reverse.mediation.calc <- function(target.gene, qtl.snp, geno, expression, annotation, covar=NULL){

  # get the expression levels of target.gene (the mediator)
  mediator.pos  <- which(colnames(expression) %in% target.gene)
  mediator.expr <- expression[, mediator.pos]
  # calculate the prob.matrix of geno
  prob.matrix <- snp.prob.matrix(geno)
  med <- NULL

  for(i in 1:ncol(expression)){
    gene.i <- colnames(expression)[i]
    target.expr <- expression[,i]
    # perform mediation analysis on gene.i with target.gene as the mediator
    med.i <- mediation.scan(target = target.expr, mediator = mediator.expr,
                            annotation = annotation[mediator.pos,], qtl.geno = prob.matrix[, qtl.snp, ],
                            covar = covar, verbose = FALSE)
    med.i <- med.i[,c("SYMBOL", "LOD")]
    med.i[1,"SYMBOL"] <- gene.i

    if(is.null(med)){
      med <- med.i
    }else{
      med <- rbind(med, med.i)
    }
  }
  med <- join(med, annotation, by = "SYMBOL")

  return(med)
}
