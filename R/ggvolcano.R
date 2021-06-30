#' Volcano plot by ggplot2
#'
#' @param data results of difference expression genes, which must contains
#'     logFC and p value.
#' @param logFC default is logFC
#' @param pvalue default if PValue
#' @param symbol default if NULL
#' @param pvalue.cutoff default is 0.01
#' @param lgf default is 2
#' @param alpha for points, default is 0.5
#' @param size for points, default is 2
#' @param colours for points, default is red, bule and gray
#' @param size.border 1
#' @param size.vhline 0.75
#' @param size.axis.text 13
#' @param size.axis.title 14
#' @param size.symbol 4
#' @param ticks.length 0.15
#' @param ticks.width 1
#' @param legend.text text for three groups, default is Up, Down and NS. NS means not significant
#' @param legend.title NULL
#' @param legend.position legend position, which can be right, left, top, bottom, or two of them. Or two numbers, eg: c(1,1)
#' @param size.legend.text 12
#' @param size.legend.title 15
#' @param symbol.pvalue.cutoff default is NULL, which takes value of pvalue.cutoff
#' @param symbol.lgf default is NULL, which takes value of lgf
#' @param cat logical, whether to output gene number for three classifications.
#' @param flip logical, default is TRUE, whether to flip coordinate
#' @param grid logical, default is FALSE, whether to show panel grid
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 aes aes_string element_line element_rect element_text geom_hline
#' @importFrom ggplot2 geom_point geom_vline ggplot scale_color_manual theme theme_bw unit
#' @importFrom ggplot2 xlab ylab coord_flip element_blank
#' @return one ggplot2 volcano picture
#' @export
#'
ggvolcano <- function(data,
                      logFC='logFC',pvalue='PValue',symbol=NULL,
                      pvalue.cutoff=0.01,lgf=2,
                      alpha=0.5,size=2,colours=c('red', "blue", "gray"),
                      size.border=1,
                      size.vhline=0.75,
                      size.axis.text=13,
                      size.axis.title=14,
                      size.symbol=4,
                      ticks.length=0.15,
                      ticks.width=1,
                      legend.text=c('Up','Down','NS'),
                      legend.title=NULL,
                      legend.position='right',
                      size.legend.text=12,
                      size.legend.title=15,
                      symbol.pvalue.cutoff=NULL,
                      symbol.lgf=NULL,
                      cat=TRUE,
                      flip=FALSE,
                      grid=FALSE){
    data[,pvalue]=-log10(data[,pvalue])
    p.cutoff <- -log10(pvalue.cutoff)
    up = data[,pvalue] >=  p.cutoff & data[,logFC] >= abs(lgf)
    down = data[,pvalue] >= p.cutoff & data[,logFC] <= -abs(lgf)
    if (length(legend.text)==1){
        legend.text=c(legend.text, "Down", "NS")
    }else if (length(legend.text)==2){
        legend.text=c(legend.text, "NS")
    }
    data$ggdirection=ifelse(up,legend.text[1],ifelse(down,legend.text[2],legend.text[3]))
    data$ggdirection=factor(data$ggdirection,levels = legend.text)
    tb=table(data$ggdirection)
    if (cat){
        cat('\n')
        cat(names(tb)[1],': ',tb[1],'\n')
        cat(names(tb)[2],': ',tb[2],'\n')
        cat(names(tb)[3],': ',tb[3],'\n')

    }
    ggdirection='only for pass cran check'
    if (length(colours)==1){
        colours=c(colours, "blue", "gray")
    }else if (length(colours)==2){
        colours=c(colours, "gray")
    }
    p <- ggplot(data = data, aes_string(x = logFC, y = pvalue)) +
        geom_point(aes(colour = ggdirection),alpha = alpha,size=size) +
        scale_color_manual(values = colours,name=legend.title)+
        theme_bw()+
        theme(legend.title = element_text(size=size.legend.title),
              legend.text = element_text(size=size.legend.text),
              panel.border = element_rect(size = size.border),
              axis.ticks.length = unit(ticks.length,'cm'),
              axis.ticks = element_line(size=ticks.width),
              axis.text = element_text(size=size.axis.text),
              axis.title = element_text(size=size.axis.title),
              legend.position = legend.position,)+
        ylab(expression(-log[10]("P Value"))) +
        xlab(expression(log[2]("Fold Change")))+
        geom_vline(xintercept = c(-abs(lgf), abs(lgf)),
                   lty = 2,
                   col = "black",
                   lwd = size.vhline)+
        geom_hline(yintercept = p.cutoff,
                   lty = 2,
                   col = "black",
                   lwd = size.vhline)
    p
    if (!is.null(symbol)){
        if (is.null(symbol.pvalue.cutoff)) symbol.pvalue.cutoff= pvalue.cutoff
        if (is.null(symbol.lgf)) symbol.lgf= lgf

        data_symbol=data[data[,pvalue] >= -log10(symbol.pvalue.cutoff) &
                             abs(data[,logFC]) >= abs(symbol.lgf), ]
        symbol='only for pass cran check'
        p<-p+ geom_text_repel(data = data_symbol,
                              aes(symbol = symbol),
                              size = size.symbol,
                              segment.color = "black")
        p
    }
    if (flip) p <- p + coord_flip()
    if (!grid) p <- p + theme(panel.grid = element_blank())
    p$direction = data$ggdirection
    p
}
