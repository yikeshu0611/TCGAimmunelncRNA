mRNA_lncRNA <- function(x){
    # mRNA/lncRNA
    message('\n使用注释文件gene code v22提取mRNA和lncRNA')
    message('网址：https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files')
    message('下载地址：https://api.gdc.cancer.gov/data/b011ee3e-14d8-4a97-aed4-e0b10f6bbe82')

    # mRNA
    message('\n\n提取 mRNA')
    cat('mRNA的基因类型是：\n              protein_coding')
    mRNA_matrix <- translation(x$expr,dic = 'mRNA')

    # lncRNA
    message('\n\n提取 lncRNA')
    type <- c("3prime_overlapping_ncRNA",
              "antisense",
              "bidirectional_promoter_lncRNA",
              "lincRNA",
              "macro_lncRNA",
              "non_coding",
              "processed_transcript",
              "sense_intronic",
              "sense_overlapping")
    cat('lncRNA的基因类型是:\n')
    cat(paste0(paste0('             ',type),collapse = '\n'))
    lncRNA_matrix <- translation(x$expr,dic = 'lncRNA')

    # group
    group <- ifelse(do::mid(rownames(mRNA_matrix),14,2) |> as.numeric() <= 9,'Tumor','Normal')
    group <- data.frame(group=group,
                        bar_code=rownames(mRNA_matrix),
                        bcr_patient_barcode=do::left(rownames(mRNA_matrix),12))

    # RNA with survival
    message('\n\n创建包含生存信息的mRNA矩阵和lncRNA矩阵')
    lncRNA_surv <- add_survial(lncRNA_matrix,x$clinical_short)
    mRNA_surv <- add_survial(mRNA_matrix,x$clinical_short)
    list(mRNA=mRNA_matrix,
         mRNA_surv=mRNA_surv,
         lncRNA=lncRNA_matrix,
         lncRNA_surv=lncRNA_surv,
         group=group,
         clinical_short=x$clinical_short
         )
}



add_survial <- function(RNA,clinical_short){
    lncnm <- colnames(RNA)
    RNA$time <- clinical_short[do::left(rownames(RNA),12),'time']
    RNA$status <- clinical_short[do::left(rownames(RNA),12),'status']
    RNA[,c('time','status',lncnm)]
}



translation <- function(expr,dic){
    if (dic=='mRNA'){
        dic_expr <- v22[v22$gene_type %in% 'protein_coding',]
    }else if (dic=='lncRNA'){
        type <- c("3prime_overlapping_ncRNA",
                  "antisense",
                  "bidirectional_promoter_lncRNA",
                  "lincRNA",
                  "macro_lncRNA",
                  "non_coding",
                  "processed_transcript",
                  "sense_intronic",
                  "sense_overlapping")
        dic_expr <- v22[v22$gene_type %in% type,]

    }

    dic_matrix <- dplyr::inner_join(x = dic_expr[,c('gene_id','gene_name')],
                                    y = expr,
                                    by='gene_id')
    dic_matrix <- dic_matrix[,colnames(dic_matrix) %not% 'gene_id']
    cat('\n共有',dic,'基因: ',nrow(dic_matrix),'个')
    if (any(duplicated(dic_matrix$gene_name))){
        dp <- dic_matrix$gene_name[duplicated(dic_matrix$gene_name)]
        cat('\n重复',length(dp),'个',', 重复的基因取中位数\n')
        pb <- txtProgressBar(max = length(unique(dp)),width = 30,style = 3)
        udp <- unique(dp)
        for (i in 1:length(udp)) {
            setTxtProgressBar(pb,i)
            row_loc <- which(dic_matrix$gene_name == udp[i])
            dp_rows <- dic_matrix[row_loc,-1]
            rp <- sapply(1:(ncol(dic_matrix)-1), function(i) median(dp_rows[,i]))
            dic_matrix[row_loc[1],2:ncol(dic_matrix)] <- rp
            dic_matrix <- dic_matrix[-(row_loc[-1]),]
        }
        cat('\n剩余',nrow(dic_matrix),'个基因')
    }
    rownames(dic_matrix) <- dic_matrix$gene_name
    dic_matrix <- dic_matrix[,colnames(dic_matrix) %not% 'gene_name']
    as.data.frame(t(dic_matrix))
}
