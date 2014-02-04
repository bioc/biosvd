##' Generates up to 5 visualizations of the eigensystem to aid in deciding which eigenfeatures and eigenassays to filter out (representing noise, steady state, steady-scale, experimental artifacts), or to aid in exploring the dynamics of expression/intensity levels over time or between different groups of assays.
##'
##' The function generates a heatmap of the eigenfeatures by assays with use of the given contrast factor (heatmap), a bar plot with the eigenexpression fractions of all eigenfeatures (fraction), a bar plot with the eigenexpression fractions of the eigenfeatures without the dominant eigenfeature(s) (zoomedFraction), the intensity levels of eigenfeatures 1 to 4 across the assays (lines), and the intensity levels of all eigenfeatures across the assays (allLines).
##' @title Generate visualizations of the eigensystem
##' @param x object of class eigensystem
##' @param contrast contrast used for the heatmap visualization (default 3)
##' @param plots plots to be shown: heatmap, fraction, zoomedFraction, lines and/or allLines
##' @param prefix string to start the plot names with (default 'biosvd')
##' @param dir directory to save the plots to (default work directory)
##' @param figure boolean indicating whether figures should be shown (TRUE) or saved as pdf (FALSE, default)
##' @return NULL
##' @import Biobase grid BiocGenerics
##' @importFrom methods callGeneric show
##' @importFrom gplots heatmap.2
##' @importFrom graphics plot
##' @export
##' @docType methods
##' @rdname plot-methods
##' @aliases plot plot,Eigensystem,ANY-method
##' @section Methods:
##' \describe{
##' \item{\code{signature(x = "Eigensystem")}}{
##' }
##' }
##' @family "eigensystem"
##' @examples
##' ## Metabolomics starvation data obtained from http://genomics-pubs.princeton.edu/StarvationMetabolomics/Download.shtml
##' data(StarvationData)
##'
##' ## Computes the eigensystem for the actual data
##' eigensystem <- compute(StarvationData)
##' ## Generates all provided plots for the eigensystem
##' plot(eigensystem)
##' ## Generates all provided plots for the eigensystem, with use of contrast 2 for the heatmap
##' plot(eigensystem, contrast=2)
##' ## Generates the fraction and lines plot for the eigensystem
##' plot(eigensystem, plots=c("fraction","lines"))
##' @references Alter O, Brown PO and Botstein D. Singular value decomposition for genome-wide expression data processing and modeling. Proc Natl Acad Sci U.S.A. 97(18), 10101-10106 (2000).
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
##'
setMethod("plot", signature(x="Eigensystem"),
          function(x, y, contrast=3,
                   plots=c("heatmap","fraction","zoomedFraction","lines","allLines"),
                   prefix="biosvd", dir=getwd(), figure=FALSE)
{
  eigensystem <- x
  plots <- match.arg(plots, several.ok=TRUE)
  .checkContrast(contrast)
  .checkPrefix(prefix)
  .checkFigure(figure)
  
  pal <- colorRampPalette(c(rgb(1,0,0), rgb(0,1,0)), space="rgb")
  nassays <- ncol(eigenassays(eigensystem))
  
  if("heatmap" %in% plots) {
    contrastMatrix <- contrast*eigenfeatures(eigensystem)
    contrastMatrix[which(contrastMatrix>1)]=1
    contrastMatrix[which(contrastMatrix<(-1))]=-1
    contrastMatrix <- (contrastMatrix - min(contrastMatrix))/(max(contrastMatrix)-min(contrastMatrix))
    if (!figure) pdf(paste(prefix,'eigenfeatureXassay.heatmap.pdf',sep="."))
    heatmap.2(contrastMatrix,
              Rowv=NA, Colv=NA, scale="none", dendrogram="none",
              col=pal, trace="none",
              xlab="Assays", ylab="Eigenfeatures",
              margins=c(9,3), main=prefix
              )
	if (!figure) dev.off()
  }

  if("fraction" %in% plots) {
    if (!figure) pdf(paste(prefix,'eigenexpressionFraction.pdf',sep="."))
    upperlimit <- ceiling(max(fractions(eigensystem))*100)/100
    barplot(rev(fractions(eigensystem)),
            width=0.8, space=0.2, horiz=TRUE,
            main=paste("Entropy = ",entropy(eigensystem),sep=""),
            xlab="Eigenexpression fraction", ylab="Eigenfeatures",
            col="red", xlim=c(0,upperlimit), cex.names=0.75, las=1
            )
    lines(c(upperlimit/4,upperlimit/4),c(0,nassays-0.5))
    lines(c(upperlimit/2,upperlimit/2),c(0,nassays-0.5))
    lines(c(upperlimit*3/4,upperlimit*3/4),c(0,nassays-0.5))
    lines(c(upperlimit,upperlimit),c(0,nassays-0.5))
    if (!figure) dev.off()
  }
  
  if("zoomedFraction" %in% plots) {
    if (!figure) pdf(paste(prefix,'eigenexpressionFractionZoom.pdf',sep="."))
    upperlimit <- ceiling(max(fractions(eigensystem)[-excludeEigenfeatures(eigensystem)])*100)/100
    barplot(rev(fractions(eigensystem)[-excludeEigenfeatures(eigensystem)]),
            width=0.8, space=0.2, horiz=TRUE,
            main=paste("Entropy = ",entropy(eigensystem),sep=""),
            xlab="Eigenexpression fraction", ylab="Eigenfeatures",
            col="red", xlim=c(0,upperlimit), cex.names=0.75, las=1
            )
    lines(c(upperlimit/4,upperlimit/4),c(0,nassays-0.5))
    lines(c(upperlimit/2,upperlimit/2),c(0,nassays-0.5))
    lines(c(upperlimit*3/4,upperlimit*3/4),c(0,nassays-0.5))
    lines(c(upperlimit,upperlimit),c(0,nassays-0.5))
    if (!figure) dev.off()
  }
  
  if("lines" %in% plots) { 
    if (!figure) pdf(paste(prefix,'eigenfeatureLevelXassays.1-4.pdf',sep="."))
    col.map <- c("black","red","blue","green")
    z <- 1
    par(mar=c(9,4,4,2))
    plot(eigenfeatures(eigensystem)[z,],
         ylim=c(-1,1), xlab="", ylab="Eigenfeature level",
         type="b", pch=16, xaxt='n', yaxt='n',
         col=col.map[z]
         )
    axis(1, at=1:ncol(eigenfeatures(eigensystem)), labels=colnames(eigenfeatures(eigensystem)), las=2)
    axis(2, at=seq(-1,1,0.5), labels=c(1,0.5,0,0.5,1), las=1)
    abline(h=0)
    for (z in c(2:4)) lines(eigenfeatures(eigensystem)[z,], pch=16, col=col.map[z], type="b")
    legend("topright", legend=c("first","second","third","fourth"), fill=col.map, bty="n")
    if (!figure) dev.off()
  }
  
  if("allLines" %in% plots) { 
    nb.plots <- floor(nassays/4)
    mod <- nassays%%4
    if (mod>1) { nb.plots <- nb.plots+1 }
    nb.rows <- ceiling(nb.plots/2)
    if (!figure) pdf(paste(prefix,'eigenfeatureLevelXassays.all.pdf',sep="."))
    par(mar=c(1,4,1,1), mfrow=c(nb.rows,2))
    for (z in c(1:nb.plots)) {
      if (z==nb.plots) {subset <- c((z*4-3):nassays)} else {subset <- c((z*4-3):(z*4))}
      col.map <- rainbow(length(subset))
      plot(eigenfeatures(eigensystem)[subset[1],],
           ylim=c(-1,1),
           xlab="", ylab="Eigenfeature level",
           type="b", pch=16, xaxt='n', yaxt='n',
           col=col.map[1]
           )
      axis(2, at=seq(-1,1,0.5), labels=c(1,0.5,0,0.5,1), las=1)
      abline(h=0)
      for (i in c(2:length(subset))) lines(eigenfeatures(eigensystem)[subset[i],], pch=16, col=col.map[i], type="b")
      legend("topright", legend=subset, fill=col.map, bty="n")
    }
    if (!figure) dev.off()
    par(mfrow=c(1,1))
  }
})