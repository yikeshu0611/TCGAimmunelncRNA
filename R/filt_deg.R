#' Filt deg by pvalue and lfc
#'
#' @param deg deg results
#' @param pvalue.cutoff 0.05 default
#' @param lfc 2 default
#'
#' @export
#'
filt_deg <- function(deg,pvalue.cutoff=0.01,lfc=2){

    if (class(deg)[1]=='DGEExact'){
        tbl <- as.data.frame(deg$table)
        p <- p.adjust(tbl$PValue, method = "BH")
        tbl$FDR <- p
        tbl$direction='NotSig'

        up = tbl[,'PValue'] <=  pvalue.cutoff & tbl[,'logFC'] >= abs(lfc)
        down = tbl[,'PValue'] <= pvalue.cutoff & tbl[,'logFC'] <= -abs(lfc)
        tbl$direction[up] <- 'Up'
        tbl$direction[down] <- 'Down'

        print(table(tbl$direction))
        tbl
    }
}
