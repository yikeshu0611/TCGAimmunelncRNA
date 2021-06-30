#' clear data by clinical_short
#'
#' @param x reslut of FPKM2df
#'
#' @export
#'
consistent <- function(x){
    for (i in 1:4) {
        x[[i]] <- x[[i]][do::left(rownames(x[[i]]),12) %in% rownames(x$clinical_short),]
    }
    x[[2]] <- cbind(x[[2]][,1:2],x[[1]])
    x[[4]] <- cbind(x[[4]][,1:2],x[[3]])
    group <- do::mid(rownames(x[[1]]),14,2) |> as.numeric()
    group <- ifelse(group <= 9,'Tumor','Control')
    group <- data.frame(group,
                        bar_code=rownames(x[[1]]),
                        bcr_patient_barcode=do::left(rownames(x[[1]]),12))
    x$group <- group
    x
}
