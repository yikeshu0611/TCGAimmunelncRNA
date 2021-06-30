library(docxR)

#   读取数据 -------

pRCC_data <- FPKM2df(trans_cart = 'F:\\immunelncRNA\\02.mRNAdownload\\gdc_download_20190703_015745.101334',
                     trans_metadata = "F:\\immunelncRNA\\02.mRNAdownload\\metadata.cart.2019-07-03.json",
                     clinical_cart = 'F:\\immunelncRNA\\05.clnical\\gdc_download_20190703_022829.964065')


#####  数据处理 ----
# 临床数据 ----
#        删除随访时间过短的病例
#
# data_filt <- filt1(d = pRCC_data,minTime = 30)
#
#
#
# age <- ifelse(pRCC_data$clinical_short$age>=60,'≥ 60','< 60')
# age <- ifelse(is.na(age),'Unknown',age)
# pRCC_data$clinical_short$age <- age
#
# #      制作临床特征三线表
# tb1 <- tableone(pRCC_data$clinical_short)
# docx_add_table(tbl = tb1,
#                file.name = '肾癌免疫相关lncRNA.docx')


# 测序数据 -------
# 提取mRNA和lncRNA
RNA_data <- mRNA_lncRNA(pRCC_data)


# 测序数据处理
# RNA_filt1 <- filt2(d = RNA_data,
#                   RNA = 'mRNA',
#                   minExpPatientRatio = 0.6,
#                   minExpmean = 1,
#                   minExpsd = 0.5)
#
#
# RNA_filt2 <- filt2(d = RNA_filt1,
#                   RNA = 'lncRNA',
#                   minExpPatientRatio = 0.6,
#                   minExpmean = 1,
#                   minExpsd = 0.5)

# 提取免疫相关的mRNA
library(msig)

M13664 <- msig_geneSymbol('immune_system_process')
immue_RNA <- immuneRNA(d = RNA_data,genes = M13664)


# correlation for lncRNA

cor <- cor_RNA(immue_RNA)
cor2 <- filt_cor(cor,minCor = 0.8,pvalue.cutoff = 0.01)

# 免疫相关lncRNA
imlncRNA <- cor2$lncRNA
imlncRNA_sample <- cor2$bar_code

# 提取数据
imlncRNA_surv <- immue_RNA$lncRNA_surv[imlncRNA_sample,
                                     c('time','status',imlncRNA)]
# 单因素

uv_res <- uv_cox(data = imlncRNA_surv,
                 time = 'time',
                 status = 'status',
                 pvalue.cutoff = 0.05)

rownames(uv_res)


# 多因素

mv_res <- mv_cox(data = imlncRNA_surv,
                 time = 'time',
                 status = 'status',
                 x = row.names(uv_res),
                 direction = 'both')
res_mulCox
final_gene <- row.names(res_mulCox)[res_mulCox$pvalue <= 0.05]

gene_riskdata <- riskdata


# 生存曲线
library(survival)
library(survminer)
fit <- survfit(Surv(time,status)~riskgroup,gene_riskdata)
ggsurvplot(fit = fit,
           data = gene_riskdata,
           pval = T,
           risk.table = T,
           legend.labs=c('High risk','Low risk'))
summary(fit,times = c(365,365*3,365*5))



# 矫正临床因素
clrs <- clinic_risk(immue_RNA,gene_riskdata)

clrs$stage <- do::Replace0(clrs$stage,c('A','B','C'))
clrs$n <- do::Replace0(clrs$n,c('a','b'))
fastStat::list.str(clrs,n = 6)

head(clrs)




uv_cox(data = clrs,time = 'time',status = 'status',drop = F)
fit <- mv_cox(data = clrs,time = 'time',
                  status = 'status')


# 多指标ROC
library(modelROC)
r <- roc(fit,x=T,method='KM')
ggplot(r,rank = TRUE)
unique(r)

# 临床相关性
clinic_index = c('stage','t','n','m')
rownames(imlncRNA_surv) <- do::left(row.names(imlncRNA_surv),12)
exprlist <- lapply(clinic_index, function(i){
    clin <- clrs[,i,drop=FALSE]
    clin <- cbind(clin,log2(imlncRNA_surv[rownames(clin),final_gene]))
    x <- reshape2::melt(data = clin,id=colnames(clin)[1])
    colnames(x)[2:3] <- c('gene','expression')
    x
})
names(exprlist) <- clinic_index

data=exprlist$stage
ggboxplot(data = data,
          x="gene",
          y="expression",
          color = colnames(data)[1],
          ylab="lncRNA expression",
          xlab='',ylim=c(-8,8)) +
    rotate_x_text(90) +
    stat_compare_means(aes_string(group=colnames(data)[1]),
                       label = "p",
                       vjust = -3)


data=exprlist$t
ggboxplot(data = data,
          x="gene",
          y="expression",
          color = colnames(data)[1],
          ylab="lncRNA expression",
          xlab='',ylim=c(-8,8)) +
    rotate_x_text(90) +
    stat_compare_means(aes_string(group=colnames(data)[1]),
                       label = "p",
                       vjust = -3)

data=exprlist$n
ggboxplot(data = data,
          x="gene",
          y="expression",
          color = colnames(data)[1],
          ylab="lncRNA expression",
          xlab='',ylim=c(-8,8)) +
    rotate_x_text(90) +
    stat_compare_means(aes_string(group=colnames(data)[1]),
                       label = "p",
                       vjust = -3)

data=exprlist$m
ggboxplot(data = data,
          x="gene",
          y="expression",
          color = colnames(data)[1],
          ylab="lncRNA expression",
          xlab='',ylim=c(-8,8)) +
    rotate_x_text(90) +
    stat_compare_means(aes_string(group=colnames(data)[1]),
                       label = "p",
                       vjust = -3)


# 主成分分析
imlncRNA_sample
imlncRNA

#  免疫相关的lncRNA
pcadata <- t(imlncRNA_surv[,-c(1,2)])

pcadata <- log2(pcadata+1)
dim(pcadata)
for (i in 1:ncol(pcadata)) {
    pcadata[,i] <- (pcadata[,i]-min(pcadata[,i]))/(max(pcadata[,i])-min(pcadata[,i]))
}
range(pcadata)

pcar <- prcomp(x = pcadata,scale = TRUE)
library(factoextra)
fviz_eig(pcar)
fviz_pca_ind(pcar,col.var = gene_riskdata$riskgroup)
fviz_pca_var(pcar,
             col.var = "contrib",repel = TRUE)


fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


predpca <- predict(pcar)
plot(x=predpca[,1],y=predpca[,2])













