#' Filt lncRNA_surv from lncRNA_surv and cor_RNA
#'
#' @param lncRNA_surv lncRNA_surv
#' @param cor results of cor_RNA, should after filt_cor()
#'
#' @export
#'
filt_lncRNA <- function(lncRNA_surv,cor){
    lncRNA_surv[cor$bar_code,c('time','status',cor$lncRNA)]
}
