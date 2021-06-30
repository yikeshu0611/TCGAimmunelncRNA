to.cph <- function (fit) {
    if (class(fit)[1] == "coxph") {
        fit <- base.rms::coxph2cph(fit)
        update(fit, x = TRUE, y = TRUE, surv = TRUE)
    }else if (class(fit)[1] == "cph") {
        update(fit, x = TRUE, y = TRUE, model = TRUE, surv = TRUE)
    }
}
