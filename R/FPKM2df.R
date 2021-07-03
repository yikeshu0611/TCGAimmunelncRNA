#' clear fpkm data and clinical data
#'
#' @param trans_cart directory of fpkm cart
#' @param trans_metadata directory of fpkm metadata
#' @param clinical_cart directory of clinical
#' @param pattern file patern defalut is FPKM.txt
#' @importFrom foreach %dopar%
#' @importFrom set %not%
#' @return mRNA, lncRNA, clinical data
#' @export
#'
FPKM2df <- function(clinical_cart,
                    trans_cart,
                    trans_metadata,
                    pattern='FPKM.txt'){

    # clinical
    cat(crayon::bgWhite('\n\n === 处理   临床数据  === '))
    cat('\nclinical是原始的临床数据，直接提取的，没有经过处理')
    clinical <- read_TCGA_Clinical(Clinical_dir = clinical_cart)


    # clinical_short
    message('\n\n提取 clinical_short')
    cat('这个是整理好的临床数据, 后面的分析都用它')
    cat(crayon::cyan(paste('\n\n共有病人: ',nrow(clinical),'个')))
    clinical_short <- clinical_short(clinical = clinical)


    cat(crayon::bgWhite('\n\n\n === 处理  测序数据    === '))
    fpkm_files <- list.files(trans_cart,pattern = pattern,full.names = T,recursive = T)
    r <- fpkm2list(fpkm_files,pattern)
    cat('\n合并各个文件成表达矩阵\n')
    expr <- fpkm_join(r)
    cat(crayon::cyan('\n样本个数: ',ncol(expr)-1))

    # trans_metadata
    fn_id <- fpkm2metadata(trans_metadata)
    colnames(expr)[-1] <- fn_id[colnames(expr)[-1],1]
    cat(crayon::cyan('\nTumor: ',sum(as.numeric(do::mid(colnames(expr)[-1],14,2)) <= 9)))
    cat(crayon::cyan('\nNormal: ',sum(as.numeric(do::mid(colnames(expr)[-1],14,2)) > 9)))
    cat(crayon::cyan('\n基因个数: ',nrow(expr)))

    list(expr=expr,
         clinical_short=clinical_short,
         clinical=clinical)
}

clinical_short <- function(clinical){
    cl <- clinical[,c("bcr_patient_barcode",
                      "days_to_death","vital_status","gender",
                      "age_at_initial_pathologic_diagnosis",
                      "pathologic_stage",
                      "pathologic_t","pathologic_n","pathologic_m")]
    cl$days_to_death <- ifelse(is.na(cl$days_to_death),
                               clinical$days_to_last_followup,
                               cl$days_to_death)
    colnames(cl) <- c("bcr_patient_barcode", "time", "status", "gender",
                      "age", "stage","t", "n", "m")
    ck <- do::NA.row.sums(cl[,c('time','status')]) >0
    cl <- cl[!ck,]
    if (anyDuplicated(cl$bcr_patient_barcode)){
        cl_check <<- cl
        stop('请联系我微信Charleszhanggo')
    }
    rownames(cl) <- cl$bcr_patient_barcode

    cl <- cl[,colnames(cl) %not% 'bcr_patient_barcode']
    cl$status <- as.numeric(ifelse(tolower(cl$status)=='dead',1,0))
    cl$time <- as.numeric(cl$time)

    cl$t <- do::Replace0(cl$t,c('a','b','c'))
    cl$t[cl$t == 'TX'] = 'Unknown'
    cl$t[is.na(cl$t)] = 'Unknown'
    cl$n[cl$n == 'NX'] = 'Unknown'
    cl$n[is.na(cl$n)] = 'Unknown'
    cl$m[cl$m == 'MX'] = 'Unknown'
    cl$m[is.na(cl$m)] = 'Unknown'
    cl$stage[is.na(cl$stage)] = 'Unknown'
    cl$age=as.numeric(cl$age)
    cat(crayon::cyan('\n收集的临床信息有: ',paste0(colnames(cl),collapse = ', ')))
    cat(crayon::cyan('删除time或者status缺失的病例：',sum(ck)))
    cat('\n生存状态status编码转换 Dead to 1, Alive to 0')
    cl
}

fpkm2metadata <- function(trans_metadata){
    json <- jsonlite::read_json(trans_metadata)
    entity_submitter_id <- sapply(1:length(json), function(i) json[[i]][["associated_entities"]][[1]][["entity_submitter_id"]])
    file_name <- sapply(1:length(json), function(i) json[[i]][["file_name"]])
    fn_id <- data.frame(entity_submitter_id)
    rownames(fn_id) <- file_name
    fn_id
}

fpkm_join <- function(r){
    pb <- txtProgressBar(max=length(r),width = 30,style = 3)
    for(i in 1:length(r)){
        setTxtProgressBar(pb = pb,value = i)
        if (i ==1){
            expr <- r[[i]]
        }else{
            expr <- dplyr::full_join(x = expr,y = r[[i]],by='gene_id')
        }
    }
    expr
}



fpkm2list <- function(fpkm_files,pattern){
    cores <- parallel::detectCores()
    cl <- snow::makeSOCKcluster(cores)
    doSNOW::registerDoSNOW(cl)
    message('\n共有 ',length(fpkm_files),' 个',pattern,' ','文件')
    pb <- txtProgressBar(max=length(fpkm_files),
                         width = 30,
                         style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    t1 <- Sys.time()
    result <-
        foreach::foreach(i=1:length(fpkm_files),
                .packages="Kendall",
                .options.snow=opts) %dopar% {
                    df <- data.table::fread(fpkm_files[i],
                                            showProgress = FALSE,
                                            data.table = FALSE)
                    colnames(df) <- c('gene_id',do::Replace0(fpkm_files[i],'.*/'))
                    df$gene_id <- gsub('\\..*','',df$gene_id)
                    df
                }
    close(pb)
    snow::stopCluster(cl)
    cat(time_segment(as.numeric(Sys.time()) - as.numeric(t1)))
    return(result)
}
