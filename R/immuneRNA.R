#' Extract immune RAN
#'
#' @param d results
#' @param genes immune genes
#'
#' @export
#'
immuneRNA <- function(d,genes){
    immGenes <- set::and(colnames(d$mRNA),genes)
    cat('共有免疫基因',length(immGenes),'个')
    d$mRNA <- d$mRNA[,immGenes]
    d$mRNA_surv <- cbind(d$mRNA_surv[,c(1,2)],d$mRNA)
    d
}
