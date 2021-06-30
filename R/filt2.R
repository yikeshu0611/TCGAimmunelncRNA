#' Filt gene expression
#'
#' @param df dfframe with subject as rows and genes as columns
#' @param minExpPatientRatio dafault is 0.4. patitent ratio of minium expression of pratient
#' @param minExpmean dafault is 0.5, minium mean gene expression
#' @param minExpsd dafault is 0.5, minium sd gene expression
#'
#' @export
#'
filt2 <- function(d,RNA='mRNA',
                 minExpPatientRatio=0.4,
                 minExpmean=0.5,
                 minExpsd=0.5){
    loc <- grepl(tolower(RNA),tolower(names(d)))
    df <- d[loc][1][[1]]
    di <- df[,colnames(df) %not% c('time','status')]
    # minExpPatientRatio
    cat('\n\n',RNA,'总数',ncol(di))

    ratio <- sapply(1:ncol(di), function(j) sum(!is.na(di[,j]))/nrow(di))
    di <- di[,ratio > minExpPatientRatio]
    cat('\n删除最少表达样本比例为 ',minExpPatientRatio,
        '的',RNA,sum(ratio < minExpPatientRatio),'个',
        '\n还剩: ', sum(ratio >= minExpPatientRatio))

    # minExpmean
    mean <- colMeans(di,na.rm = T)
    di <- di[,mean > minExpmean]
    cat('\n\n删除表达量过低的,如均值小于',minExpmean,'的',RNA,
        sum(mean < minExpmean),'个')
    cat('\n还剩: ',sum(mean >= minExpmean))

    # minExpsd
    sd <- sapply(1:ncol(di), function(j) sd(di[,j],T))
    di <- di[,sd > minExpsd]
    cat('\n\n删除表达过于一致的,如标准差小于',minExpsd,'的',RNA,
        sum(sd < minExpsd),'个')
    cat('\n还剩: ',sum(sd >= minExpsd))

    d[loc][1][[1]] <- di
    d[loc][2][[1]] <- cbind(d[loc][2][[1]][,c(1,2)],di)
    d
}
