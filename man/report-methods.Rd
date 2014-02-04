\docType{methods}
\name{report,Eigensystem-method}
\alias{report}
\alias{report,Eigensystem-method}
\title{Generates a summarizing report of the eigensystem}
\usage{
\S4method{report}{Eigensystem}(x, eigenfeature.xaxis = 2,
  eigenfeature.yaxis = 1, colorIdAssays = rep(1, ncol(matrix(eigensystem))),
  colorIdFeatures = rep(1, nrow(matrix(eigensystem))), contrast = 3,
  prefix = "biosvd", dir = getwd())
}
\arguments{
  \item{x}{output of compute, plot or exclude}

  \item{eigenfeature.xaxis}{positive integer giving the
  eigenfeature to be shown on the x-axis of the polar
  plots}

  \item{eigenfeature.yaxis}{positive integer giving the
  eigenfeature to be shown on the y-axis of the polar
  plots}

  \item{colorIdAssays}{vector of integers, factor or a
  string variable in the pheno data of the original data
  set in case from class ExpressionSet, with assay
  annotation used for coloring purposes on the polar plot}

  \item{colorIdFeatures}{vector of integers, factor or a
  string variable in the feature data of the original data
  set in case from class ExpressionSet, with feature
  annotation used for coloring purposes on the polar plot}

  \item{contrast}{positive integer giving the contrast used
  for the heatmap visualization}

  \item{prefix}{prefix to start the plot and html file
  with}

  \item{dir}{directory to which the plots and html table
  are saved}
}
\description{
Generates a summarizing report of the eigensystem
}
\details{
The function generates an html report of the eigensystem,
containing polar plots for the assays and features
separately that show the assays/features according to their
correlation with two eigenfeatures/eigenassays, and a table
with the list of features, sortable according to their
coordinates, radius and phase in the polar plot.

The function also generates a visualization of the sorted
data according to the two given eigenfeatures and
eigenassays, showing a heatmap of the features by assays
with colored feature/assay annotation information when
provided, a heatmap of the features by eigenassays, and the
intensity levels of the two sorted eigenassays across all
features.
}
\section{Methods}{
  \describe{ \item{\code{signature(x = "Eigensystem")}}{ }
  }
}
\examples{
## Metabolomics starvation data obtained from http://genomics-pubs.princeton.edu/StarvationMetabolomics/Download.shtml
data(StarvationData)

## Computes the eigensystem for the actual data
eigensystem <- compute(StarvationData)
## Exclude the eigenfeatures representing steady-state intensity
eigensystem <- exclude(eigensystem)
## Computes the eigensystem on the variance in the data after filtering out stead-state intensity
eigensystem <- compute(eigensystem, apply="variance")
## No exclusion of eigenfeatures representing steady-scale variance
eigensystem <- exclude(eigensystem, excludeEigenfeatures=0)

## Generate report for eigenfeatures 1 and 2 without coloring of the assays and features
report(eigensystem)
## Generate report for eigenfeatures 2 and 3
report(eigensystem, eigenfeature.xaxis=2, eigenfeature.yaxis=3)
## Generate report for eigenfeatures 1 and 2 with coloring of the assays as indicated by the variable "Starvation" in the ExpressionSet {\\tt StarvationData} (starvation: carbon vs. nitrogen)
report(eigensystem, colorIdAssays="Starvation")
## Generate report for eigenfeatures 1 and 2 with use of prefix "AwesomeStudy" for the names of the plots and html table
report(eigensystem, prefix="AwesomeStudy")
}
\author{
Anneleen Daemen \email{daemen.anneleen@gene.com}, Matthew
Brauer \email{brauer.matthew@gene.com}
}
\seealso{
Other X.eigensystem.: \code{\link{compute}},
\code{\link{compute,data.frame-method}},
\code{\link{compute,Eigensystem-method}},
\code{\link{compute,Eigensystem-method}},
\code{\link{compute,ExpressionSet-method}},
\code{\link{compute,matrix-method}}; \code{\link{exclude}},
\code{\link{exclude,Eigensystem-method}},
\code{\link{exclude,Eigensystem-method}};
\code{\link{plot}},
\code{\link{plot,Eigensystem,ANY-method}},
\code{\link{plot,Eigensystem,ANY-method}}
}

