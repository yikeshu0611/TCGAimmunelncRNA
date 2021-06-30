#' filt cor_RNA results
#'
#' @param cor results of cor_RNA()
#' @param minCor min correlation, default is 0.4
#' @param pvalue.cutoff pvalue cutoff ,default is 0.05
#'
#' @export
#'
filt_cor <- function(cor,minCor=0.4,pvalue.cutoff=0.05){
    cr <- cor$cor
    cr <- suppressMessages(do::complete.data(cr))
    ck <- abs(cr$cor) >= minCor &  cr$pvalue<= pvalue.cutoff
    x <- cr[ck,]
    cor$cor <- x
    cor$mRNA <- unique(x[,1])
    cor$lncRNA <- unique(x[,2])

    cat('mRNA: ',length(unique(x[,1])))
    cat('\nlncRNA: ',length(unique(x[,2])))
    cor
}
