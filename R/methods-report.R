##' Generates a summarizing report of the eigensystem
##'
##' The function generates an html report of the eigensystem, containing polar plots for the assays and features separately that show the assays/features according to their correlation with two eigenfeatures/eigenassays, and a table with the list of features, sortable according to their coordinates, radius and phase in the polar plot.
##' 
##' The function also generates a visualization of the sorted data according to the two given eigenfeatures and eigenassays, showing a heatmap of the features by assays with colored feature/assay annotation information when provided, a heatmap of the features by eigenassays, and the intensity levels of the two sorted eigenassays across all features.
##' @title Generates a summarizing report of the eigensystem
##' @param x output of compute, plot or exclude
##' @param eigenfeature.xaxis positive integer giving the eigenfeature to be shown on the x-axis of the polar plots
##' @param eigenfeature.yaxis positive integer giving the eigenfeature to be shown on the y-axis of the polar plots
##' @param colorIdAssays vector of integers, factor or a string variable in the pheno data of the original data set in case from class ExpressionSet, with assay annotation used for coloring purposes on the polar plot
##' @param colorIdFeatures vector of integers, factor or a string variable in the feature data of the original data set in case from class ExpressionSet, with feature annotation used for coloring purposes on the polar plot
##' @param contrast positive integer giving the contrast used for the heatmap visualization
##' @param prefix prefix to start the plot and html file with
##' @param dir directory to which the plots and html table are saved
##' @return NULL
##' @import Biobase grid BiocGenerics
##' @importFrom methods callGeneric show
##' @importFrom gplots heatmap.2
##' @importFrom hwriter hwrite hwriteImage
##' @importFrom ReportingTools HTMLReport publish finish
##' @export
##' @docType methods
##' @rdname report-methods
##' @aliases report report,Eigensystem-method
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
##' ## Exclude the eigenfeatures representing steady-state intensity
##' eigensystem <- exclude(eigensystem)
##' ## Computes the eigensystem on the variance in the data after filtering out stead-state intensity
##' eigensystem <- compute(eigensystem, apply="variance")
##' ## No exclusion of eigenfeatures representing steady-scale variance
##' eigensystem <- exclude(eigensystem, excludeEigenfeatures=0)
##'
##' ## Generate report for eigenfeatures 1 and 2 without coloring of the assays and features
##' report(eigensystem)
##' ## Generate report for eigenfeatures 2 and 3
##' report(eigensystem, eigenfeature.xaxis=2, eigenfeature.yaxis=3)
##' ## Generate report for eigenfeatures 1 and 2 with coloring of the assays as indicated by the variable "Starvation" in the ExpressionSet {\tt StarvationData} (starvation: carbon vs. nitrogen)
##' report(eigensystem, colorIdAssays="Starvation")
##' ## Generate report for eigenfeatures 1 and 2 with use of prefix "AwesomeStudy" for the names of the plots and html table
##' report(eigensystem, prefix="AwesomeStudy") 
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
setMethod("report", "Eigensystem", 
          function(x, eigenfeature.xaxis=2, eigenfeature.yaxis=1,
                   colorIdAssays=rep(1,ncol(matrix(eigensystem))),
                   colorIdFeatures=rep(1,nrow(matrix(eigensystem))),
                   contrast=3, prefix="biosvd", dir=getwd())
{
  eigensystem <- x
  .checkEigenfeatureAxes(eigenfeature.xaxis, eigenfeature.yaxis, eigensystem)
  colorIdAssays <- .checkColorID(colorIdAssays, eigensystem, "assays")
  colorIdFeatures <- .checkColorID(colorIdFeatures, eigensystem, "features")
  .checkContrast(contrast)
  .checkPrefix(prefix)
  
  eigensystem.sorted <- sort(eigensystem, decreasing=FALSE, eigenfeature.xaxis, eigenfeature.yaxis, colorIdFeatures)
  nr.matrix <- nrow(matrix(eigensystem))
  
  # Heatmap of sorted features by assay
  pal <- colorRampPalette(c(rgb(1,0,0), rgb(0,1,0)), space="rgb")
  contrastMatrix <- contrast*matrix(eigensystem.sorted)
  contrastMatrix[which(contrastMatrix>1)] <- 1
  contrastMatrix[which(contrastMatrix<(-1))] <- -1
  contrastMatrix <- (contrastMatrix - min(contrastMatrix))/(max(contrastMatrix)-min(contrastMatrix))
  unique.col.ids <- sort(unique(colorIdAssays), na.last=NA)
  nb.col.ids <- length(unique.col.ids)
  col.assays <- rep(0,ncol(matrix(eigensystem.sorted)))
  col.map <- rainbow(nb.col.ids)
  for (z in c(1:nb.col.ids)) col.assays[which(colorIdAssays %in% unique.col.ids[z])] <- col.map[z]
  unique.row.ids <- sort(unique(colorIdFeatures(eigensystem.sorted)), na.last=NA)
  nb.row.ids <- length(unique.row.ids)
  col.features <- rep(0,nrow(matrix(eigensystem.sorted)))
  row.map <- rainbow(nb.row.ids)
  for (z in c(1:nb.row.ids)) {col.features[which(colorIdFeatures(eigensystem.sorted) %in% unique.row.ids[z])] <- row.map[z]}
  pdf(paste(prefix,'sortedfeatureXassay.heatmap.pdf',sep="."))
  heatmap.2(contrastMatrix,
            Rowv=NA, Colv=NA,
            RowSideColors=col.features, ColSideColors=col.assays,
            scale="none", dendrogram="none", col=pal, trace="none",
            xlab="Assays", ylab="Features",
            labRow=NA, margins=c(9,3),
            main=prefix, key=TRUE
            )
  legend("left",
         legend=c("Assay annotation",as.character(unique.col.ids),"","Feature annotation",as.character(unique.row.ids)),
         fill=c("white",col.map,"white","white",row.map),
         bty="n", border=FALSE, cex=0.7, y.intersp=0.7
         )
  dev.off()
  
  # Heatmap of sorted features by eigenassays
  contrast <- contrast*10
  contrastMatrix <- contrast*eigenassays(eigensystem.sorted)
  contrastMatrix[which(contrastMatrix>1)]=1
  contrastMatrix[which(contrastMatrix<(-1))]=-1
  contrastMatrix <- (contrastMatrix - min(contrastMatrix))/(max(contrastMatrix)-min(contrastMatrix))
  pdf(paste(prefix,'sortedfeatureXeigenassay.heatmap.pdf',sep="."))
  heatmap.2(contrastMatrix,
            Rowv=NA, Colv=NA,
            scale="none", dendrogram="none",
            col=pal, tracecol=rgb(0.5,0.5,0), trace="none",
            xlab="Eigenassays", ylab="Features",
            labRow=NA, margins=c(3,3), main=prefix
            )
  dev.off()

  # Sorted features by levels
  pdf(paste(prefix,'sortedfeatureXlevel.pdf',sep="."))
  plot(eigenassays(eigensystem.sorted)[,1],c(nr.matrix:1),
       xlab="Feature level", ylab="Features",
       type="l", col="red",
       xlim=c(floor(100*min(eigenassays(eigensystem.sorted)[,1:2]))/100, ceiling(100*max(eigenassays(eigensystem.sorted)[,1:2]))/100)
       )
  abline(v=0)
  lines(eigenassays(eigensystem.sorted)[,2],c(nr.matrix:1), col="green", type="l")
  dev.off()
  
  # Polar plot of assays in function of 2 eigenfeatures
  coordinates.assays <- base::matrix(0,nrow=ncol(matrix(eigensystem)),ncol=2)
  for (z in c(1:ncol(eigenassays(eigensystem)))) {
    scale.factor <- sqrt(matrix(eigensystem)[,z] %*% matrix(eigensystem)[,z])
    coordinates.assays[z,] <- c(assaycorrelations(eigensystem)[eigenfeature.xaxis,z]/scale.factor, assaycorrelations(eigensystem)[eigenfeature.yaxis,z]/scale.factor)
  }
  radii.assays <- signif(sqrt(coordinates.assays[,1]^2+coordinates.assays[,2]^2),3)
  names(radii.assays) <- colnames(matrix(eigensystem))
  phase.assays <- atan(assaycorrelations(eigensystem)[eigenfeature.yaxis,]/assaycorrelations(eigensystem)[eigenfeature.xaxis,])/pi
  names(phase.assays) <- colnames(matrix(eigensystem))
  coordinates.assays <- signif(coordinates.assays,3)
  rownames(coordinates.assays) <- colnames(matrix(eigensystem))
  
  pdf(paste(prefix,'polarplot.assays.eigenfeature',eigenfeature.xaxis,'vs',eigenfeature.yaxis,'pdf',sep="."))
  vp0 <- viewport(x=0,width=0.05,just="left",name="vp0")
  vp1 <- viewport(x=0.1,y=0.1,width=0.75,height=0.75,just=c("left","bottom"),name="vp1")
  vp2 <- viewport(x=0.1,y=0,width=0.75,height=0.1,just=c("left","bottom"),name="vp2")
  vp3 <- viewport(x=1,width=0.2,just="right",name="vp3")
  pushViewport(vp0)
  grid.text(paste("Assay correlation with eigenassay ",eigenfeature.yaxis,sep=""), y=0.5, rot=90)
  upViewport()
  pushViewport(vp1)
  grid.circle(x=0.5,y=0.5,r=0.5,gp=gpar(lty="dashed"))
  grid.circle(x=0.5,y=0.5,r=0.25,gp=gpar(lty="dashed",fill="grey"))
  grid.lines(x=unit(c(0,1),"npc"),y=unit(c(0.5,0.5),"npc"),arrow=NULL)
  grid.lines(x=unit(c(0.5,0.5),"npc"),y=unit(c(0,1),"npc"),arrow=NULL)
  for (z in c(1:nb.col.ids)) {
    indices <- which(colorIdAssays %in% unique.col.ids[z])
    grid.points(x=unit((coordinates.assays[indices,1]+1)/2,"npc"),y=unit((coordinates.assays[indices,2]+1)/2,"npc"), pch=z,gp=gpar(col=col.map[z],cex=0.7))
  }
  grid.text(c(1:length(colorIdAssays)),just="left",x=(coordinates.assays[,1]+1)/2+0.02,y=(coordinates.assays[,2]+1)/2)
  grid.lines(x=unit(c(0.5,(coordinates.assays[1,1]+1)/2),"npc"),y=unit(c(0.5,(coordinates.assays[1,2]+1)/2),"npc"),arrow=arrow(angle=30, length=unit(0.02,"npc"),ends="last",type="open"))
  upViewport()
  pushViewport(vp2)
  grid.text(paste("Assay correlation with eigenassay ",eigenfeature.xaxis,sep=""), x=0.5, y=0.5)
  upViewport()
  pushViewport(vp3)
  grid.points(pch=1:nb.col.ids,x=unit(rep(0.5,nb.col.ids),"lines"),y=unit(1,"npc")-unit(c(1:nb.col.ids),"lines"),gp=gpar(col=col.map))
  grid.text(unique.col.ids,just="left",x=unit(rep(1.5,nb.col.ids),"lines"), y=unit(1,"npc")-unit(c(1:nb.col.ids),"lines"),gp=gpar(col=col.map))
  upViewport()
  dev.off()
  
  # Polar plot of features in function of 2 eigenfeatures
  coordinates.features <- base::matrix(0,nrow=nr.matrix,ncol=2)
  for (z in c(1:nr.matrix)) {
    scale.factor <- sqrt(matrix(eigensystem)[z,] %*% matrix(eigensystem)[z,])
    coordinates.features[z,] <- c(featurecorrelations(eigensystem)[eigenfeature.xaxis,z]/scale.factor, featurecorrelations(eigensystem)[eigenfeature.yaxis,z]/scale.factor)
  }
  radii.features <- signif(sqrt(coordinates.features[,1]^2+coordinates.features[,2]^2),3)
  names(radii.features) <- rownames(matrix(eigensystem))
  phase.features <- atan(featurecorrelations(eigensystem)[eigenfeature.yaxis,]/featurecorrelations(eigensystem)[eigenfeature.xaxis,])/pi
  names(phase.features) <- rownames(matrix(eigensystem))
  coordinates.features <- signif(coordinates.features,3)
  rownames(coordinates.features) <- rownames(matrix(eigensystem))
  phase.features.converted <- phase.features*pi
  phase.features.converted[which(coordinates.features[,1]<0)] <- phase.features.converted[which(coordinates.features[,1]<0)]+pi
  phase.features.converted[which(coordinates.features[,1]>0 && coordinates.features[,2]<0)] <- phase.features.converted[which(coordinates.features[,1]>0 && coordinates.features[,2]<0)]+(2*pi)
  phase.features.converted <- signif(phase.features.converted,3)

  pdf(paste(prefix,'polarplot.features.eigenfeature',eigenfeature.xaxis,'vs',eigenfeature.yaxis,'pdf',sep="."))
  vp0 <- viewport(x=0,width=0.05,just="left",name="vp0")
  vp1 <- viewport(x=0.1,y=0.1,width=0.75,height=0.75,just=c("left","bottom"),name="vp1")
  vp2 <- viewport(x=0.1,y=0,width=0.75,height=0.1,just=c("left","bottom"),name="vp2")
  vp3 <- viewport(x=1,width=0.2,just="right",name="vp3")
  pushViewport(vp0)
  grid.text(paste("Feature correlation with eigenfeature ",eigenfeature.yaxis,sep=""), y=0.5, rot=90)
  upViewport()
  pushViewport(vp1)
  grid.circle(x=0.5,y=0.5,r=0.5,gp=gpar(lty="dashed"))
  grid.circle(x=0.5,y=0.5,r=0.25,gp=gpar(lty="dashed",fill="grey"))
  grid.lines(x=unit(c(0,1),"npc"),y=unit(c(0.5,0.5),"npc"),arrow=NULL)
  grid.lines(x=unit(c(0.5,0.5),"npc"),y=unit(c(0,1),"npc"),arrow=NULL)
  for (z in c(1:nb.row.ids)) {
    indices <- which(colorIdFeatures %in% unique.row.ids[z])
    grid.points(x=unit((coordinates.features[indices,1]+1)/2,"npc"),y=unit((coordinates.features[indices,2]+1)/2,"npc"), pch=z,gp=gpar(col=row.map[z],cex=0.7))
  }
  grid.lines(x=unit(c(0.5,(coordinates.features[1,1]+1)/2),"npc"),y=unit(c(0.5,(coordinates.features[1,2]+1)/2),"npc"),arrow=arrow(angle=30, length=unit(0.02,"npc"),ends="last",type="open"))
  upViewport()
  pushViewport(vp2)
  grid.text(paste("Feature correlation with eigenfeature ",eigenfeature.xaxis,sep=""), x=0.5, y=0.5)
  upViewport()
  pushViewport(vp3)
  grid.points(pch=1:nb.row.ids,x=unit(rep(0.5,nb.row.ids),"lines"),y=unit(1,"npc")-unit(c(1:nb.row.ids),"lines"),gp=gpar(col=row.map))
  grid.text(unique.row.ids,just="left",x=unit(rep(1.5,nb.row.ids),"lines"), y=unit(1,"npc")-unit(c(1:nb.row.ids),"lines"),gp=gpar(col=row.map))
  upViewport()
  dev.off()

  # Html summary report with polarplots
  polarplot.info <- cbind(rownames(matrix(eigensystem)),radii.features,phase.features.converted,coordinates.features[,c(2,1)])
  colnames(polarplot.info)=c("Feature","Radius","Phase",paste("Coordinate to eigenfeature ",eigenfeature.yaxis,sep=""),paste("Coordinate to eigenfeature ",eigenfeature.xaxis,sep=""))
  if (class(featureMatrix(eigensystem))!="NULL" && ncol(featureMatrix(eigensystem)) > 0) polarplot.info <- cbind(polarplot.info, featureMatrix(eigensystem))
  
  myPage <- HTMLReport(baseUrl=dir,
                       reportDirectory=prefix,
                       shortName=paste("report.eigenfeature",eigenfeature.xaxis,"vs",eigenfeature.yaxis,"html",sep="."),
                       title="SVD analysis"
                       )
  publish(hwrite(paste("SVD analysis: ",prefix,sep=""), heading=2), myPage)
  polarplotsmat <- rbind(c("Assay polar plot","Feature polar plot"), hwriteImage(c(paste(prefix,'polarplot.assays.eigenfeature',eigenfeature.xaxis,'vs',eigenfeature.yaxis,'pdf',sep="."), paste(prefix,'polarplot.features.eigenfeature',eigenfeature.xaxis,'vs',eigenfeature.yaxis,'pdf',sep=".")), table=FALSE))
  publish(hwrite(polarplotsmat, heading=3, br=TRUE, center=TRUE), myPage)
  publish(hwrite(''), myPage)
  publish(hwrite(paste("Feature list for ",prefix,sep=""), heading=3), myPage)
  publish(as.data.frame(polarplot.info), myPage)
  finish(myPage)
})