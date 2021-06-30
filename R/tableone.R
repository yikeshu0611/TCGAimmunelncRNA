#' create table one for classified variable
#'
#' @param data survival data
#'
#' @export
#'
tableone <- function(data){
    ck <- sapply(1:ncol(data), function(i) is.character(data[,i]) | is.factor(data[,i]))
    d1 <- data[,ck,drop=FALSE]
    nrow <- sum(sapply(1:ncol(d1), function(i) length(unique(d1[,i]))))


    for (i in 1:ncol(d1)) {
        if (i==1){
            j=0
            matrix(rep('',4*nrow),ncol=4,dimnames = list(NULL,
                                                         c('Characteristics','Variable',
                                                           'Entire dataset','Percentages (%)'))) |>
                as.data.frame() -> df
        }
        udi <- unique(d1[,i])
        udi <- sort(udi)
        for (k in 1:length(udi)) {
            j = j +1
            if (k==1) df[j,1] <- colnames(d1)[i]

            df[j,2] <- udi[k]
            df[j,3] <- sum(d1[,i] == udi[k])
            df[j,4] <- round(sum(d1[,i] == udi[k])/length(d1[,i])*100,2)
        }
    }
    df
}
