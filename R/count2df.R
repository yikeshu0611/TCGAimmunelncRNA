#' clear count data and clinical data
#'
#' @param count_cart directory of count cart
#' @param count_metadata directory of count metadata
#' @param clinical_cart directory of clinical
#' @param pattern file patern defalut is htseq.counts.gz
#' @return mRNA, lncRNA, clinical data
#' @export
count2df <- function(count_cart,
                     count_metadata,
                     clinical_cart,
                     pattern='htseq.counts.gz'){
    FPKM2df(count_cart,
            count_metadata,
            clinical_cart,
            pattern)
}
