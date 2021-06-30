year_survRate <- function(survfit,times){
    x <- summary(fit,times = times)
    data.frame(model = deparse(substitute(survfit)),
               strata = x$strata,
               times = x$time,
               surv = x$surv,
               lower = x$lower,
               upper = x$upper,
               n.risk = x$n.risk)
}
