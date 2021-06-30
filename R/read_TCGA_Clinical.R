#' Read Clinical Data of TCGA
#'
#' @param Clinical_dir directory of TCGA clinical data, which contains one or more xml files
#' @importFrom magrittr %>%
#' @return one dataframe
#' @export
#'
read_TCGA_Clinical <- function(Clinical_dir){
    xml_files <- list.files(path = Clinical_dir,
                            pattern = 'xml',
                            recursive = TRUE,
                            full.names = TRUE)
    cat('\n提取xml文件个数: ',length(xml_files),'\n')
    pb <- txtProgressBar(min = 0,max = length(xml_files),initial = 0,
                         width = 30,style = 3)
    xmls <- lapply(1:length(xml_files), function(i){
        setTxtProgressBar(pb = pb,value = i)
        xml_i(link = xml_files[i])
    })
    close(pb)
    r <- do.call(plyr::rbind.fill, xmls)
    r
}


xml_i <- function(link){
    children <- xml2::read_html(link) %>%
        do::all_children()

    nodes_name <- children %>%
        xml2::xml_name()
    nodes_text <- children %>%
        xml2::xml_text()

    text_zero <- nchar(nodes_text) > 0
    nodes_text <- nodes_text[text_zero]
    nodes_name <- nodes_name[text_zero]

    lastest <- !rev(duplicated(rev(nodes_name)))
    nodes_name <- nodes_name[lastest]
    nodes_text <- nodes_text[lastest]

    mt <- matrix(nodes_text,nrow = 1,dimnames = list(NULL, nodes_name))
    data.frame(mt,check.names = FALSE)
}
