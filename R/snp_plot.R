#' SNP frequency plot
#'
#'@description Function to draw a SNP frequency plot. Requires ggplot2.
#'
#' @param x A VCF dataframe that contains the SNP data from a single isolate only.
#'   vcf_transformation.
#' @param col Column of the dataframe which is used to colorize each SNP position.
#' @param genome.length The length of the reference genome. If set to default the last SNP position is used.
#' @param mark.repeats Dataframe that contains the repeat regions of the genome.
#' @param mark.lowGQ Numeric value >= 0
#' @param mark.lowQUAL Numeric value >= 0
#' @return A list of two elements. First element is the input dataframe with SNP
#'   specificities added. The second element is a counting talbe with all SNP
#'   specificity combinations.
#' @examples
#' iso1 <- bac[bac$ISO == "iIsolate1",]
#' bacsnp.plot(iso1, col = "SPEC")
#' bacsnp.plot(iso1, col = "SPEC", genome.length = 125000)
#' @import ggplot2
#' @export
bacsnp.plot <- function(x, col = "SPEC",
                        genome.length = NULL,
                        mark.repeats = NULL,
                        mark.lowGQ = NULL,
                        mark.lowQUAL = NULL){

  x <- x[!is.na(x$REL.ALT1),]

  if(is.numeric(genome.length)){
    b <- seq(0, genome.length, by = 5000)
    b[1] <- 1
    b <- append(b, genome.length)
    #b[which(b == max(b))] <- genome.length
  }else{
    genome.length <- max(x$POS)
    b <- seq(0, genome.length, by = 5000)
    b[1] <- 1
    b <- append(b, genome.length)
  }


  #if(mark.repeats == TRUE){
  if(is.data.frame(mark.repeats) == TRUE){
    include_repeats <- geom_rect(data = mark.repeats, aes(xmin = mark.repeats$xstart, xmax = mark.repeats$xend, ymin = -Inf, ymax = Inf), alpha = 0.8, fill = "grey")
  }else{
    include_repeats <- geom_rect()
  }


  if(is.numeric(mark.lowGQ)){
    lowGQ <- x[x$GQ < mark.lowGQ,]
    x <- x[x$GQ >= mark.lowGQ,]
    mark_GQ <- geom_point(data = lowGQ, mapping = aes(x = as.numeric(lowGQ$POS), y = lowGQ$REL.ALT1), shape = 0, color = "lightgrey")
    #mark_GQ <- geom_point(data = lowGQ, colour = 'grey', shape = 8)
  }else{
    mark_GQ <- geom_point()
  }

  if(is.numeric(mark.lowQUAL)){
    lowQUAL <- x[x$QUAL < mark.lowQUAL,]
    x <- x[x$QUAL >= mark.lowQUAL,]
    mark_QUAL <- geom_point(data = lowQUAL, mapping = aes(x = as.numeric(lowQUAL$POS), y = lowQUAL$REL.ALT1), shape = 2, color = "lightgrey")
    #mark_QUAL <- geom_point(data = lowQUAL, colour = 'black', shape = 8)
  }else{
    mark_QUAL <- geom_point()
  }


  p <- ggplot() +

    theme_bw() +
    scale_x_continuous(expand = c(0.01, 0), breaks = b) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(colour = "black", fill= NA, size = 1)
          #legend.position = "right",
          #legend.title = element_text("SNP specificity")
    ) +
    include_repeats +
    mark_GQ +
   # theme(legend.position = "bottom") +
    mark_QUAL +
    geom_point(data = x, mapping = aes(x = as.numeric(x$POS), y = x$REL.ALT1, color = x[,col])) +
    ggtitle(as.character(unique(x$ISO))) +
    xlab("position") +
    ylab("relativ frequency")

  return(p)

}

