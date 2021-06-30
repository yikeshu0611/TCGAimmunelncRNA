#' Filt subject by minium time
#'
#' @param d result of FPKM2df()
#' @param minTime default is 30
#'
#' @export
#'
filt1 <- function(d,minTime=30){
    bar0 <- rownames(d$clinical_short)
    cat(crayon::cyan('\n总病例数: ',nrow(d$clinical_short)))
    ck <- !is.na(d$clinical_short$time) & !is.na(d$clinical_short$status)
    if (sum(!ck)>0){
        cat(crayon::cyan('\n\n删除随访时间time或者生存状态status为空值的病例 ',sum(!ck),' 个病例'))
        cat('\n还剩 ',sum(ck),' 个病例')
        d$clinical_short <- d$clinical_short[ck,]
    }
    ck <- rownames(d$clinical_short) %in% do::left(colnames(d$expr)[-1],12)
    if (sum(!ck) >0){
        cat(crayon::cyan('\n\n删除没有进行测序的病例: ',sum(!ck),' 个病例'))
        cat('\n还剩: ',sum(ck))
        d$clinical_short <- d$clinical_short[ck,]
    }
    # minExpPatientRatio
    d$clinical_short$time <- as.numeric(d$clinical_short$time)

    ck <- d$clinical_short$time >= minTime
    if (sum(!ck)>0){
        cat(crayon::cyan('\n删除随访时间小于',minTime,'的病例: ',sum(!ck),' 个病例'))
        cat('\n病例最终剩余 ',sum(ck),' 个病例')
        d$clinical_short <- d$clinical_short[ck,]
    }

    ck <- do::left(colnames(d$expr)[-1],12) %in% rownames(d$clinical_short)
    if (sum(!ck)>0){
        cat('\n\n\n总样本: ',length(ck),' 个')
        cat(crayon::cyan('\n相应地, 删除测序数据中没有临床信息的样本',sum(!ck),'个'))
        cat('\n剩余 ',sum(ck),' 个样本')
        d$expr <- d$expr[,c(TRUE,ck)]

        cla <- ifelse(do::mid(colnames(d$expr)[-1],14,2) |> as.numeric() <= 9,'Tumor','Normal')
        cat('\nTumor: ',sum(cla=='Tumor'))
        cat('\nNormal: ',sum(cla=='Normal'))
    }
    d
}
