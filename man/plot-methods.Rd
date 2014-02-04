\docType{methods}
\name{plot,Eigensystem,ANY-method}
\alias{plot}
\alias{plot,Eigensystem,ANY-method}
\title{Generate visualizations of the eigensystem}
\usage{
\S4method{plot}{Eigensystem,ANY}(x, y, contrast = 3, plots = c("heatmap",
  "fraction", "zoomedFraction", "lines", "allLines"), prefix = "biosvd",
  dir = getwd(), figure = FALSE)
}
\arguments{
  \item{x}{object of class eigensystem}

  \item{contrast}{contrast used for the heatmap
  visualization (default 3)}

  \item{plots}{plots to be shown: heatmap, fraction,
  zoomedFraction, lines and/or allLines}

  \item{prefix}{string to start the plot names with
  (default 'biosvd')}

  \item{dir}{directory to save the plots to (default work
  directory)}

  \item{figure}{boolean indicating whether figures should
  be shown (TRUE) or saved as pdf (FALSE, default)}

  \item{y}{the y coordinates of points in the plot,
  \emph{optional} if \code{x} is an appropriate structure.}
}
\description{
Generates up to 5 visualizations of the eigensystem to aid
in deciding which eigenfeatures and eigenassays to filter
out (representing noise, steady state, steady-scale,
experimental artifacts), or to aid in exploring the
dynamics of expression/intensity levels over time or
between different groups of assays.
}
\details{
The function generates a heatmap of the eigenfeatures by
assays with use of the given contrast factor (heatmap), a
bar plot with the eigenexpression fractions of all
eigenfeatures (fraction), a bar plot with the
eigenexpression fractions of the eigenfeatures without the
dominant eigenfeature(s) (zoomedFraction), the intensity
levels of eigenfeatures 1 to 4 across the assays (lines),
and the intensity levels of all eigenfeatures across the
assays (allLines).
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
## Generates all provided plots for the eigensystem
plot(eigensystem)
## Generates all provided plots for the eigensystem, with use of contrast 2 for the heatmap
plot(eigensystem, contrast=2)
## Generates the fraction and lines plot for the eigensystem
plot(eigensystem, plots=c("fraction","lines"))
}
\author{
Anneleen Daemen \email{daemen.anneleen@gene.com}, Matthew
Brauer \email{brauer.matthew@gene.com}
}
\references{
Alter O, Brown PO and Botstein D. Singular value
decomposition for genome-wide expression data processing
and modeling. Proc Natl Acad Sci U.S.A. 97(18), 10101-10106
(2000).
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
\code{\link{report}},
\code{\link{report,Eigensystem-method}},
\code{\link{report,Eigensystem-method}}
}

