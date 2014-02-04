\docType{methods}
\name{sort,Eigensystem-method}
\alias{sort}
\alias{sort,Eigensystem-method}
\title{Sorts the data according to two eigenfeatures and eigenassays}
\usage{
\S4method{sort}{Eigensystem}(x, decreasing = FALSE, eigenfeature.xaxis = 2,
  eigenfeature.yaxis = 1, colorIdFeatures = rep(1,
  nrow(matrix(eigensystem))))
}
\arguments{
  \item{x}{object of class eigensystem}

  \item{decreasing}{FALSE}

  \item{eigenfeature.xaxis}{first eigenfeature to sort
  against (default 2)}

  \item{eigenfeature.yaxis}{second eigenfeature to sort
  against (default 1)}

  \item{colorIdFeatures}{vector of integers or a string
  variable in the feature data of the original data set in
  case from class ExpressionSet, with feature annotation
  used for coloring purposes on the polar plot}
}
\value{
Object of class eigensystem
}
\description{
Sorts the data according to two eigenfeatures and
eigenassays
}
\details{
The function sorts the data according to two certain
eigenfeatures and eigenassays. This gives a global picture
of the dynamics of expression/intensities, in which
individual features and assays are classified in groups of
similar regulation and function or similar cellular state
and biological phenotype.
}
\section{Methods}{
  \describe{ \item{\code{signature(x = "Eigensystem")}}{ }
  }
}
\examples{
## Data obtained from http://genomics-pubs.princeton.edu/StarvationMetabolomics/Download.shtm
data(StarvationData)

## Computes the eigensystem for the actual data
eigensystem <- compute(StarvationData)
## Excludes the eigenfeatures representing steady-state intensity
eigensystem <- exclude(eigensystem)

## Sort the data according to eigenfeature 1 and 2
eigensystem.sorted <- sort(eigensystem)
## Sort the data according to eigenfeature 3 and 4
eigensystem.sorted <- sort(eigensystem, eigenfeature.xaxis=4, eigenfeature.yaxis=3)
## Visualization of the data after sorting to eigenfeature 3 and 4
plot(eigensystem.sorted)
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

