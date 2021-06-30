#' scale for dataframe
#' (x - mean(x))/sd(x)
#' @param data dataframe
#'
#' @export
#'
fpkm_scale <- function(data){
    v <- colnames(data) %not% c('time','status')
    for (i in v) {
        data[,i]=(data[,i]-min(data[,i]))/(max(data[,i])-min(data[,i]))
    }
    data
}
