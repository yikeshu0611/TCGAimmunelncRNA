#' Looping for Unix Cox Regression
#'
#' @param data data
#' @param time time variable
#' @param status status variable
#' @param x variable names for univariable cox regression. If missing, it will be column names of data except y and adjust
#' @param adjust adjust variable names for univariable cox regression
#' @param round digital round, 3 is defaulted
#' @param pvalue.cutoff threshold for p value to show star. 0.05 is defaulted
#' @param drop logical, whether to drop data by p_value
#'
#' @return univariable cox regression results
#' @export
uv_cox <- function (data, time, status, x, adjust,
                    round = 3,
                    pvalue.cutoff = 0.05,drop=TRUE){
    if (is.factor(data[, status])) stop(status, " shoud be numeric not factor")
    if (missing(adjust)) {
        if (missing(x)) x = colnames(data) %not% c(time, status)
        pb <- txtProgressBar(max = length(x),width = 30,style = 3)
        for (i in 1:length(x)) {
            setTxtProgressBar(pb,i)
            if (i == 1) {
                res = NULL
                class = NULL
            }
            if (is.factor(data[, x[i]])) {
                rep_len = length(levels(data[, x[i]])) -
                    1
                if (rep_len > 1) {
                    rep_class = rep(i + 1, rep_len)
                    class = c(class, rep_class)
                }
                else {
                    class = c(class, 0)
                }
            }else {
                class = c(class, 0)
            }
            vari <- x[i]
            if (make.names(vari) != vari) vari = paste0('`',vari,'`')
            formu = paste0("survival::Surv(", time, ",", status,
                           ")~", vari)
            cox.i = survival::coxph(as.formula(formu), data = data)
            cox.sum = summary(cox.i)
            cox_coef1 = as.data.frame(cox.sum$coefficients)[,c(3, 4, 5)]
            cox_coef2 = as.data.frame(cox.sum$conf.int)[, c(1,3, 4)]
            res.cbind = cbind(cox_coef2, cox_coef1)
            res.i = round(res.cbind, round)
            res = rbind(res, res.i)
        }
    }else {
        if (missing(x))
            x = colnames(data) %not% c(time, status, adjust)
        pb <- txtProgressBar(max = length(x),width = 30,style = 3)
        for (i in 1:length(x)) {
            setTxtProgressBar(pb,i)
            if (i == 1) {
                res = NULL
                class = NULL
            }
            if (is.factor(data[, x[i]])) {
                rep_len = length(levels(data[, x[i]])) -
                    1
                if (rep_len > 1) {
                    rep_class = rep(i + 1, rep_len)
                    class = c(class, rep_class)
                }
                else {
                    class = c(class, 0)
                }
            }
            else {
                class = c(class, 0)
            }
            if (i == 1) {
                formu = paste0("survival::Surv(", time, ",",
                               status, ")~", paste0(adjust, collapse = "+"))
                cox.i = survival::coxph(as.formula(formu), data = data)
                cox.sum = summary(cox.i)
                cox_coef1 = as.data.frame(cox.sum$coefficients)
                nub_row = nrow(cox_coef1)
            }
            vari <- x[i]
            if (make.names(vari) != vari) vari = paste0('`',vari,'`')
            formu = paste0("survival::Surv(", time, ",", status,
                           ")~", paste0(c(adjust, vari), collapse = "+"))
            cox.i = survival::coxph(as.formula(formu), data = data)
            cox.sum = summary(cox.i)
            cox_coef1 = as.data.frame(cox.sum$coefficients)[-(1:nub_row),
                                                            c(3, 4, 5)]
            cox_coef2 = as.data.frame(cox.sum$conf.int)[-(1:nub_row),
                                                        c(1, 3, 4)]
            res.cbind = cbind(cox_coef2, cox_coef1)
            res.i = round(res.cbind, round)
            res = rbind(res, res.i)
        }
    }
    res$risk <- ifelse(res$`exp(coef)` >1,'High','Low')
    if (drop) res <- res[res$`Pr(>|z|)` <= pvalue.cutoff,]
    message('\n',sum(res$`Pr(>|z|)` <= pvalue.cutoff),' significant genes')
    colnames(res)=c('OR','Low95','High95','se','z','pvalue','risk')
    rownames(res) <- gsub('`','',rownames(res))
    res
}
