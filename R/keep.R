#' keep values
#'
#' @param ... one or more value
#'
#' @export
#'
keep <- function(...){
    keep <- do::get_names(...)
    leave <- set::not(ls(envir = .GlobalEnv),keep)

    parse(text = sprintf('rm(%s,envir = .GlobalEnv)',
            paste0(leave,collapse = ','))) |>
        eval()

}
