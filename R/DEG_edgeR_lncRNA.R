#' DEG analysis using edgeR package for mRNA
#'
#' @param d results of FPKM2df()
#'
#' @export
#'
DEG_edgeR_lncRNA <- function(d){
    dl <- edgeR::DGEList(counts = t(d$lncRNA),
                         group = d$group$group)
    message('Normalization')
    dl <- edgeR::calcNormFactors(dl)
    message('Estimate Common Negative Binomial Dispersion ')
    dl <- edgeR::estimateCommonDisp(dl)
    message('Estimate Empirical Bayes Tagwise Dispersion Values')
    dl <- edgeR::estimateTagwiseDisp(dl)
    message('Exact Tests for Differences between Two Groups ')
    edgeR::exactTest(dl,pair = unique(d$group$group))
}
