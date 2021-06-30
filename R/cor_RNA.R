#' correlation of mRNA and lncRNA
#'
#' @param data data of FPKM2df()
#' @importFrom foreach %dopar%
#' @export
#'
cor_RNA <- function(data){
    mRNA <- data$mRNA[data$group$bar_code[data$group=='Tumor'],]
    message('\nuse ',nrow(mRNA), ' tumors for correlation analysis')
    message('mRNA:',ncol(mRNA))
    lncRNA <- data$lncRNA[data$group$bar_code[data$group=='Tumor'],]
    message('lncRNA:',ncol(lncRNA))
    cores <- parallel::detectCores()
    cl <- snow::makeSOCKcluster(cores)
    doSNOW::registerDoSNOW(cl)
    pb <- utils::txtProgressBar(max=ncol(mRNA),
                         width = 30,
                         style=3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    result <-
        foreach::foreach(i=1:ncol(mRNA),
                .packages="Kendall",
                .options.snow=opts) %dopar% {
                    x <- as.data.frame(t(sapply(1:ncol(lncRNA), function(j) as.data.frame(t(do.call('c',cor.test(mRNA[,i],lncRNA[,j]))))[,c(1,4,3)])))
                    colnames(x) <- c('t','cor','pvalue')
                    x <- do::numeric.it(x,c('t','cor','pvalue'))
                    x <- round(x,3)
                    x$mRNA <- colnames(mRNA)[i]
                    x$lncRNA <- colnames(lncRNA)
                    x
                }
    close(pb)
    snow::stopCluster(cl)
    co <- do.call(rbind,result)
    co$regulate <- ifelse(co$cor>0,'positive','negative')
    cor <- co[,c("mRNA", "lncRNA", "t", "cor", "pvalue",'regulate')]
    list(cor=cor,
         bar_code=rownames(mRNA),
         mRNA=colnames(mRNA),
         lncRNA=colnames(lncRNA))
}


