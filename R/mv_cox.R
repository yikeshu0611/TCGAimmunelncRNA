#' multi-variable for cox regression
#'
#' @param data data
#' @param time time variable
#' @param status status variable
#' @param x x variable
#' @param direction should be one of no, both, backward, forward
#' @param ... other arg
#'
#' @importFrom survival coxph Surv
#' @export
#'
mv_cox <- function (data, time, status, x, direction = "no",
          ...){
    if (is.factor(data[, status])) stop(status, " shoud be numeric not factor")
    if (missing(x)) x = colnames(data) %not% c(time, status)
    xdata <- data[,x,drop=FALSE]
    string <- 'coxph(Surv(%s,%s)~%s,data=%s)'
    for (i in 1:length(x)) {
        if (x[i] != make.names(x[i])){
            xi=x[i]
            x[i] = sprintf('`%s`',x[i])
            # colnames(data)[colnames(data) == xi] <- x[i]
        }
    }
    x
    colnames(data)
    s2 <- sprintf(string,time,status,
                   paste0(x, collapse = "+"),
                   as.character(substitute(data)))
    cox <- eval(parse(text = s2))

    if (direction != "no"){
        cox = step(object = cox, direction = direction)
    }
    df1 <- as.data.frame(matrix(round(summary(cox)[["conf.int"]][,c(1,3,4)],3),ncol = 3))
    df2 <- as.data.frame(matrix(round(summary(cox)[["coefficients"]][,c(3,4,5)],3),ncol = 3))
    x <- data.frame(cbind(df1,df2))
    row.names(x)=do::Replace0(row.names(summary(cox)[["conf.int"]]),'`')
    colnames(x)=c("OR", "Low95", "High95", "se", "z", "pvalue")
    x
    riskscore=round(predict(cox,newdata = data,type = 'risk'),3)
    riskgroup=ifelse(riskscore>median(riskscore),'High','Low')

    riskdata=data.frame(riskscore=riskscore,
                        riskgroup=riskgroup,
                        time=data[,time],
                        status=data[,status])
    riskdata=cbind(riskdata,xdata)
    row.names(riskdata)=row.names(data)
    res_mulCox <<- x
    riskdata <<- riskdata
    cox
}
