#' time segment for seconds
#'
#' @param x seconds
#'
#' @return character
#' @export
#'
time_segment <- function(x){
    p=ifelse(x<1*60,paste0(round(x,2),' seconds'),
             ifelse(x<1*60*60,paste0(round(x/60,2),' minutes'),
                    ifelse(x<1*60*60*24,paste0(round(x/60/60,2),' hours'),
                           ifelse(x<1*60*60*24*365,paste0(round(x/60/60/24,2),' days'),
                                  ifelse(x<1*60*60*24*365*100,paste0(round(x/60/60/24/365,2),' years'),
                                         paste0(round(x/60/60/24/365/100,2),' century'))))))
    p
}
