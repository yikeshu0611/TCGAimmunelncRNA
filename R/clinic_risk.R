#' cbind clinical data and riskdata
#'
#' @param immue_RNA immue_RNA list data
#' @param riskdata riskdata dataframe
#'
#' @export
#'
clinic_risk <- function(immue_RNA,riskdata){
    clinical <- immue_RNA$clinical_short
    clinical <- clinical[rowSums(clinical == 'Unknown',na.rm = T) == 0,]
    ck <- complete.cases(clinical[do::left(rownames(riskdata),12),])
    d1 <- clinical[do::left(rownames(riskdata),12),][ck,]
    d1$sample <- do::left(rownames(d1),12)
    rd <- riskdata[,'riskscore',drop=FALSE]
    rd$sample <- do::left(rownames(riskdata),12)
    m <- merge(x = d1,y = rd,by='sample')
    rownames(m) <- m$sample
    m[,colnames(m) %not% 'sample']
}
