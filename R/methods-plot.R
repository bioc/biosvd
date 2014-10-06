##' Generates up to ten visualizations of the eigensystem to aid in deciding which eigenfeatures and eigenassays to filter out (representing noise, steady state, steady-scale, experimental artifacts), or to aid in exploring the dynamics of expression/intensity levels over time or between different groups of assays.
##'
##' The function generates a heatmap of the eigenfeatures by assays with use of the given contrast factor (eigenfeatureHeatmap), a heatmap of the features by eigenassays with use of the given contrast factor (eigenassayHeatmap), a heatmap of the features by assays, with features sorted according to two selected eigenfeatures (sortedHeatmap), a bar plot with the eigenexpression fractions of all eigenfeatures (fraction), a screeplot for the eigenexpression fractions (scree), a bar plot with the eigenexpression fractions of the eigenfeatures without the dominant eigenfeature(s) (zoomedFraction), the intensity levels of selected eigenfeatures across the assays (by default eigenfeatures 1-4) (lines), the intensity levels of all eigenfeatures across the assays (allLines), polar plot for the assays according to their correlation with two eigenassays (eigenassayPolar), and polar plot for the features according to their correlation with two eigenfeatures (eigenfeaturePolar).
##' @title Generate visualizations of the eigensystem
##' @param x object of class Eigensystem
##' @param y object of class EigensystemPlotParam
##' @param ... Additional arguments that can be passed on
##' @return NULL
##' @import Biobase grid BiocGenerics
##' @importFrom methods callGeneric show
##' @importFrom graphics plot
##' @importFrom NMF aheatmap
##' @export
##' @docType methods
##' @rdname plot-methods
##' @aliases plot plot,Eigensystem,EigensystemPlotParam-method
##' @section Methods:
##' \describe{
##' \item{\code{signature(x = "Eigensystem", y = "EigensystemPlotParam")}}{
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
##' params <- new("EigensystemPlotParam")
##' if (.Platform$OS.type != "windows") plot(eigensystem, params)
##' ## Generates all provided plots for the eigensystem, with use of contrast 2 for the heatmap
##' contrast(params) <- 2
##' if (.Platform$OS.type != "windows") plot(eigensystem, params)
##' ## Generates the fraction and lines plot for the eigensystem
##' params <- new("EigensystemPlotParam")
##' plots(params) <- c("fraction","lines")
##' if (.Platform$OS.type != "windows") plot(eigensystem, params)
##' @references Alter O, Brown PO and Botstein D. Singular value decomposition for genome-wide expression data processing and modeling. Proc Natl Acad Sci U.S.A. 97(18), 10101-10106 (2000).
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
##'
setMethod("plot", signature(x="Eigensystem", y="EigensystemPlotParam"),
          function(x, y, ...) {
            eigensystem <- x
            params <- y
            plots <- match.arg(plots(params), c("eigenfeatureHeatmap","eigenassayHeatmap","sortedHeatmap","fraction","scree","zoomedFraction","lines","allLines","eigenassayPolar","eigenfeaturePolar"), several.ok=TRUE)
            contrast <- .checkContrast(contrast(params))
            prefix <- .checkPrefix(prefix(params))
            figure <- .checkFigure(figure(params))
            dims <- .checkDimensions(eigensystem, params)
            polarAxes <- .checkEigenfeatureAxes(eigensystem, whichPolarAxes(params))
            colormaps <- .checkColorMaps(eigensystem, params)
            assayColorMap(params) <- colormaps$assayColorMap
            featureColorMap(params) <- colormaps$featureColorMap
            annotation.colors <- unlist(colormaps, recursive=FALSE)
            names(annotation.colors) <- gsub("assayColorMap.","",names(annotation.colors))
            names(annotation.colors) <- gsub("featureColorMap.","",names(annotation.colors))
            annAssays <- as.data.frame(assayMatrix(eigensystem)[,names(colormaps$assayColorMap)[length(colormaps$assayColorMap) > 0]])
            names(annAssays) <- names(colormaps$assayColorMap)[length(colormaps$assayColorMap) > 0]
            annFeatures <- as.data.frame(featureMatrix(eigensystem)[,names(colormaps$featureColorMap)[length(colormaps$featureColorMap) > 0]])
            names(annFeatures) <- names(colormaps$featureColorMap)[length(colormaps$featureColorMap) > 0]
            pal <- palette(params)
            prefix <- paste0(path(params), "/", prefix)
         
            ########## eigenfeature heatmap ##########
            if("eigenfeatureHeatmap" %in% plots) {
              contrastMatrix <- contrast * eigenfeatures(eigensystem)
              contrastMatrix[contrastMatrix > 1] <- 1
              contrastMatrix[contrastMatrix < -1] <- -1

              if(negativeValues(params)) {
                contrastMatrix <- contrastMatrix/max(contrastMatrix)
              } else {
                contrastMatrix <- (contrastMatrix - min(contrastMatrix))/(max(contrastMatrix)-min(contrastMatrix))
              }
              
              if (!figure) pdf(paste(prefix,'eigenfeatureXassay.heatmap.pdf',sep="."))
              aheatmap(contrastMatrix,
                       Rowv=NA, Colv=NA,
                       annCol=annAssays,
                       annColors=annotation.colors,
                       color=pal(100))
              if (!figure) dev.off()
            }

            ########## eigenassay heatmap ##########
            if("eigenassayHeatmap" %in% plots) {
              contrastMatrix <- contrast * eigenassays(eigensystem)
              contrastMatrix[contrastMatrix > 1] <- 1
              contrastMatrix[contrastMatrix < -1] <- -1

              if(negativeValues(params)) {
                contrastMatrix <- contrastMatrix/max(contrastMatrix)
              } else {
                contrastMatrix <- (contrastMatrix - min(contrastMatrix))/(max(contrastMatrix)-min(contrastMatrix))
              }
              
              if (!figure) pdf(paste(prefix,'featureXeigenassay.heatmap.pdf',sep="."))
              aheatmap(contrastMatrix,
                       Rowv=NA, Colv=NA,
                       annRow=annFeatures,
                       annColors=annotation.colors,
                       color=pal(100))
              if (!figure) dev.off()
            }

            ########## sorted heatmap ##########
            if("sortedHeatmap" %in% plots) {
              contrastMatrix <- contrast * matrix(eigensystem)
              contrastMatrix[contrastMatrix > 1] <- 1
              contrastMatrix[contrastMatrix < -1] <- -1

              if(negativeValues(params)) {
                contrastMatrix <- contrastMatrix/max(contrastMatrix)
              } else {
                contrastMatrix <- (contrastMatrix - min(contrastMatrix))/(max(contrastMatrix)-min(contrastMatrix))
              }
              
              coordinates <- project(eigensystem, axes=polarAxes, type="feature")
              polar.order <- order(coordinates[,"theta"], decreasing=FALSE)
              contrastMatrix <- contrastMatrix[polar.order,]
              annFeatures <- annFeatures[polar.order, , drop=FALSE]
              
              if (!figure) pdf(paste(prefix,'featureXassay.sortedheatmap.pdf',sep="."))
              aheatmap(contrastMatrix,
                       Rowv=NA, Colv=NA,
                       annRow=annFeatures,
                       annCol=annAssays,
                       annColors=annotation.colors,
                       color=pal(100))
              if (!figure) dev.off()
            }

            ########## scree ##########
            if("scree" %in% plots) {
              if (!figure) pdf(paste(prefix,'eigenexpressionFraction.pdf',sep="."))
              upperlimit <- ceiling(max(fractions(eigensystem))*100)/100
              plot(fractions(eigensystem),
                      main=paste("Entropy = ", entropy(eigensystem),sep=""),
                      ylab="Eigenexpression fraction", xlab="Eigenfeatures",
                      pch=16, col="red", ylim=c(0,upperlimit)
                      )
              if (!figure) dev.off()
            }

            ########## fraction ##########
            if("fraction" %in% plots) {
              if (!figure) pdf(paste(prefix,'eigenexpressionFraction.pdf',sep="."))
              upperlimit <- ceiling(max(fractions(eigensystem))*100)/100
              barplot(rev(fractions(eigensystem)),
                      width=0.8, space=0.2, horiz=TRUE,
                      main=paste("Entropy = ", entropy(eigensystem),sep=""),
                      xlab="Eigenexpression fraction", ylab="Eigenfeatures",
                      col="red", xlim=c(0,upperlimit), cex.names=0.75, las=1
                      )
              lines(c(upperlimit/4,upperlimit/4),c(0,dims$nassays-0.5))
              lines(c(upperlimit/2,upperlimit/2),c(0,dims$nassays-0.5))
              lines(c(upperlimit*3/4,upperlimit*3/4),c(0,dims$nassays-0.5))
              lines(c(upperlimit,upperlimit),c(0,dims$nassays-0.5))
              if (!figure) dev.off()
            }

            ########## zoomed fraction ##########
            if("zoomedFraction" %in% plots) {
              if (!figure) pdf(paste(prefix,'eigenexpressionFractionZoom.pdf',sep="."))
              upperlimit <- ceiling(max(fractions(eigensystem)[-excludeEigenfeatures(eigensystem)])*100)/100
              barplot(rev(fractions(eigensystem)[-excludeEigenfeatures(eigensystem)]),
                      width=0.8, space=0.2, horiz=TRUE,
                      main=paste("Entropy = ", entropy(eigensystem),sep=""),
                      xlab="Eigenexpression fraction", ylab="Eigenfeatures",
                      col="red", xlim=c(0,upperlimit), cex.names=0.75, las=1
                      )
              lines(c(upperlimit/4,upperlimit/4),c(0,dims$nassays-0.5))
              lines(c(upperlimit/2,upperlimit/2),c(0,dims$nassays-0.5))
              lines(c(upperlimit*3/4,upperlimit*3/4),c(0,dims$nassays-0.5))
              lines(c(upperlimit,upperlimit),c(0,dims$nassays-0.5))
              if (!figure) dev.off()
            }
            
            ########## line plot, selected eigenfeatures ##########
            if("lines" %in% plots) { 
              if (!figure) pdf(paste(prefix,'eigenfeatureLevelXassays.selected.pdf',sep="."))

              to.be.plotted <- eigenfeatures(eigensystem)[whichEigenfeatures(params),]
              col.map <- c("#000000FF",  rainbow(length(whichEigenfeatures(params))-1));

              par(mar=c(9,4,4,2))
              plot(to.be.plotted[1,],
                   ylim=c(-1,1), xlab="", ylab="Eigenfeature level",
                   type="n", pch=16, xaxt='n', yaxt='n', lwd=4,
                   col=col.map[1]
                   )
              axis(1, at=1:ncol(to.be.plotted), labels=colnames(to.be.plotted), las=2)
              axis(2, at=seq(-1,1,0.5), labels=c(1,0.5,0,0.5,1), las=1)
              abline(h=0)
              for (z in length(whichEigenfeatures(params)):1) {
                lines(to.be.plotted[z,], pch=16, col=col.map[z], type="b", lwd=4)
              }
              legend("topright", legend= whichEigenfeatures(params), fill=col.map, bty="n")
              if (!figure) dev.off()
            }

            ########## line plot, all eigenfeatures ##########            
            if("allLines" %in% plots) {
              
              if (!figure) pdf(paste(prefix,'eigenfeatureLevelXassays.all.pdf',sep="."))
              oldpar <- par()
              par(mar=c(1,4,1,1), mfrow=c(ceiling(ceiling(dims$nassays/4)/2), 2))
              
              lines.per.plot <- 4
              col.map <- c("#000000FF",  rainbow(lines.per.plot-1));

              plotnum <- rep(1:ceiling(dims$nassays/lines.per.plot), each=lines.per.plot)[1:dims$nassays]
              plots <- tapply(1:dims$nassays, plotnum,
                              function(ef) {
                                plot(eigenfeatures(eigensystem)[ef[1],],
                                     ylim=c(-1,1),
                                     xlab="", ylab="Eigenfeature level",
                                     type="n", pch=16, xaxt='n', yaxt='n',
                                     col=col.map[1]
                                     )
                                axis(2, at=seq(-1,1,0.5), labels=c(1,0.5,0,0.5,1), las=1)
                                abline(h=0)
                                for (i in rev(ef)) {
                                  lines(eigenfeatures(eigensystem)[i,], pch=16, col=col.map[((i-1)  %% lines.per.plot) + 1], type="b")
                                }
                                legend("topright", legend=ef, fill=col.map, bty="n")
                              })
              par(mar=oldpar$mar, mfrow=oldpar$mfrow)
              if (!figure) dev.off()
            }

            ########## eigenAssay polar plot ##########
            if("eigenassayPolar" %in% plots) {

              coordinates <- project(eigensystem, polarAxes, type="assays")
              
              if (!figure) pdf(paste(prefix,'eigenassayPolar',polarAxes[1],'vs',polarAxes[2],'pdf',sep="."))
              
              vp0 <- viewport(x=0,width=0.05,just="left",name="vp0")
              vp1 <- viewport(x=0.1,y=0.1,width=0.75,height=0.75,just=c("left","bottom"),name="vp1")
              vp2 <- viewport(x=0.1,y=0,width=0.75,height=0.1,just=c("left","bottom"),name="vp2")
              vp3 <- viewport(x=1,width=0.2,just="right",name="vp3")
              pushViewport(vp0)
              grid.text(paste("Assay correlation with eigenassay ",polarAxes[2],sep=""), y=0.5, rot=90)
              upViewport()
              pushViewport(vp1)
              grid.circle(x=0.5,y=0.5,r=0.5,gp=gpar(lty="dashed"))
              grid.circle(x=0.5,y=0.5,r=0.25,gp=gpar(lty="dashed",fill="grey"))
              grid.lines(x=unit(c(0,1),"npc"),y=unit(c(0.5,0.5),"npc"),arrow=NULL)
              grid.lines(x=unit(c(0.5,0.5),"npc"),y=unit(c(0,1),"npc"),arrow=NULL)

              if (length(colormaps$assayColorMap)==0) {
              	grid.points(x=unit((coordinates[,1]+1)/2,"npc"),y=unit((coordinates[,2]+1)/2,"npc"), pch=18,gp=gpar(col="black",cex=0.7))
              } else {
                for (z in 1:length(colormaps$assayColorMap[[1]])) {
                  indices <- which(annAssays[,1]==names(colormaps$assayColorMap[[1]])[z])
                  grid.points(x=unit((coordinates[indices,1]+1)/2,"npc"),y=unit((coordinates[indices,2]+1)/2,"npc"), pch=z,gp=gpar(col=colormaps$assayColorMap[[1]][z],cex=0.7))
                }              	
              }
              grid.text(c(1:nrow(coordinates)),just="left",x=(coordinates[,1]+1)/2+0.02,y=(coordinates[,2]+1)/2)
              grid.lines(x=unit(c(0.5,(coordinates[1,1]+1)/2),"npc"),y=unit(c(0.5,(coordinates[1,2]+1)/2),"npc"),arrow=arrow(angle=30, length=unit(0.02,"npc"),ends="last",type="open"))
              upViewport()
              pushViewport(vp2)
              grid.text(paste("Assay correlation with eigenfeature ",polarAxes[1],sep=""), x=0.5, y=0.5)
              upViewport()
              pushViewport(vp3)
              if (length(colormaps$assayColorMap)>0) {
                grid.points(pch=1:length(colormaps$assayColorMap[[1]]),x=unit(rep(0.5,length(colormaps$assayColorMap[[1]])),"lines"),y=unit(1,"npc")-unit(c(1:length(colormaps$assayColorMap[[1]])),"lines"),gp=gpar(col=colormaps$assayColorMap[[1]]))
                grid.text(names(colormaps$assayColorMap[[1]]),just="left",x=unit(rep(1.5,length(colormaps$assayColorMap[[1]])),"lines"), y=unit(1,"npc")-unit(c(1:length(colormaps$assayColorMap[[1]])),"lines"),gp=gpar(col=colormaps$assayColorMap[[1]]))
    	      }
              upViewport()
              
              if (!figure) dev.off()
            }
            
            ########## eigenFeature polar plot ##########
            if("eigenfeaturePolar" %in% plots) {

              coordinates <- project(eigensystem, polarAxes, type="features")
              
              if (!figure) pdf(paste(prefix,'eigenfeaturePolar',polarAxes[1],'vs',polarAxes[2],'pdf',sep="."))
              
              vp0 <- viewport(x=0,width=0.05,just="left",name="vp0")
              vp1 <- viewport(x=0.1,y=0.1,width=0.75,height=0.75,just=c("left","bottom"),name="vp1")
              vp2 <- viewport(x=0.1,y=0,width=0.75,height=0.1,just=c("left","bottom"),name="vp2")
              vp3 <- viewport(x=1,width=0.2,just="right",name="vp3")
              pushViewport(vp0)
              grid.text(paste("Feature correlation with eigenfeature ",polarAxes[2],sep=""), y=0.5, rot=90)
              upViewport()
              pushViewport(vp1)
              grid.circle(x=0.5,y=0.5,r=0.5,gp=gpar(lty="dashed"))
              grid.circle(x=0.5,y=0.5,r=0.25,gp=gpar(lty="dashed",fill="grey"))
              grid.lines(x=unit(c(0,1),"npc"),y=unit(c(0.5,0.5),"npc"),arrow=NULL)
              grid.lines(x=unit(c(0.5,0.5),"npc"),y=unit(c(0,1),"npc"),arrow=NULL)

              if (length(colormaps$featureColorMap)==0) {
              	grid.points(x=unit((coordinates[,1]+1)/2,"npc"),y=unit((coordinates[,2]+1)/2,"npc"), pch=18,gp=gpar(col="black",cex=0.7))
              } else {
                for (z in 1:length(colormaps$featureColorMap[[1]])) {
                  indices <- which(annFeatures[,1]==names(colormaps$featureColorMap[[1]])[z])
                  grid.points(x=unit((coordinates[indices,1]+1)/2,"npc"),y=unit((coordinates[indices,2]+1)/2,"npc"), pch=z,gp=gpar(col=colormaps$featureColorMap[[1]][z],cex=0.7))
                }              	
              }
              grid.lines(x=unit(c(0.5,(coordinates[1,1]+1)/2),"npc"),y=unit(c(0.5,(coordinates[1,2]+1)/2),"npc"),arrow=arrow(angle=30, length=unit(0.02,"npc"),ends="last",type="open"))
              upViewport()
              pushViewport(vp2)
              grid.text(paste("Feature correlation with eigenfeature ",polarAxes[1],sep=""), x=0.5, y=0.5)
              upViewport()
              pushViewport(vp3)
              if (length(colormaps$featureColorMap)>0) {
                grid.points(pch=1:length(colormaps$featureColorMap[[1]]),x=unit(rep(0.5,length(colormaps$featureColorMap[[1]])),"lines"),y=unit(1,"npc")-unit(c(1:length(colormaps$featureColorMap[[1]])),"lines"),gp=gpar(col=colormaps$featureColorMap[[1]]))
                grid.text(names(colormaps$featureColorMap[[1]]),just="left",x=unit(rep(1.5,length(colormaps$featureColorMap[[1]])),"lines"), y=unit(1,"npc")-unit(c(1:length(colormaps$featureColorMap[[1]])),"lines"),gp=gpar(col=colormaps$featureColorMap[[1]]))
    	      }
              upViewport()
              
              if (!figure) dev.off()
            }

          })
