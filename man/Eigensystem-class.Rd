\docType{class}
\name{Eigensystem-class}
\alias{Eigensystem}
\alias{Eigensystem-class}
\alias{apply}
\alias{assayMatrix}
\alias{assaycorrelations}
\alias{colorIdFeatures}
\alias{eigenassays}
\alias{eigenexpressions}
\alias{eigenfeatures}
\alias{entropy}
\alias{excludeEigenfeatures}
\alias{featureMatrix}
\alias{featurecorrelations}
\alias{fractions}
\alias{matrix}
\alias{show,Eigensystem-method}
\alias{signMatrix}
\title{The Eigensystem class}
\description{
\code{Eigensystem} is a list-based class for storing the
results of applying Singular Value Decomposition (SVD) to a
feature by assay data set. Objects are normally created by
\code{\link{compute,Eigensystem-method}}.
}
\details{
The \code{Eigensystem} class stores the original data and
all SVD-derived information obtained with \code{compute}.
Data in the \code{Eigensystem} are organized into different
slots, \sQuote{matrix}, \sQuote{signMatrix},
\sQuote{assayMatrix}, \sQuote{featureMatrix},
\sQuote{eigenassays}, \sQuote{eigenexpressions},
\sQuote{eigenfeatures}, \sQuote{assaycorrelations},
\sQuote{featurecorrelations}, \sQuote{fractions},
\sQuote{entropy}, \sQuote{apply},
\sQuote{excludeEigenfeatures}, and
\sQuote{colorIdFeatures}. Brief descriptions of these slots
are provided below.
}
\section{Slots}{
  \code{Eigensystem} objects contain the following slots
  \describe{ \item{\code{matrix}:}{\code{matrix} containing
  the feature by assay data without missing values}
  \item{\code{signMatrix}:}{\code{matrix} containing the
  sign of each element in matrix}
  \item{\code{assayMatrix}:}{\code{matrix} containing
  additional information about assays, with rows as assays
  and columns as additional variables}
  \item{\code{featureMatrix}:}{\code{matrix} containing
  additional information about features, with rows as
  features and columns as additional variables}
  \item{\code{eigenassays}:}{\code{matrix} containing the
  feature by eigenassay data, with each column in
  eigenassays corresponding to a left singular vector,
  representing genome-wide expression, proteome-wide
  abundance or metabolome-wide intensity in the
  corresponding eigenassay}
  \item{\code{eigenexpressions}:}{\code{numeric} vector
  containing the eigenexpression fraction of each
  {eigenfeature, eigenassay}-pair, constituting the
  diagonal elements of the diagonal matrix connecting the
  left and right singular values; the diagonal matrix
  reflects the decoupling and decorrelation of the data,
  with expression of each eigenfeature restricted to the
  corresponding eigenassay}
  \item{\code{eigenfeatures}:}{\code{matrix} containing the
  eigenfeatures by assay data, with each row corresponding
  to a right singular vector, representing the expression,
  abundance or intensity of the corresponding eigenfeature
  across all assays}
  \item{\code{assaycorrelations}:}{\code{matrix} containing
  the correlation between the eigenassays as rows and the
  assays as columns}
  \item{\code{featurecorrelations}:}{\code{matrix}
  containing the correlation between the eigenfeatures as
  rows and features as columns}
  \item{\code{fractions}:}{\code{numeric} vector containing
  the eigenexpression fraction for each {eigenfeature,
  eigenassay}-pair, defined as the relative fraction of
  overall expression that each eigenfeature and eigenassay
  capture} \item{\code{entropy}:}{\code{numeric} value
  between 0 and 1 giving the Shannon entropy as measure for
  data complexity, with an entropy of 0 corresponding to an
  ordered and redundant data set with all expression
  captured by a single {eigenfeature, eigenassay}-pair, and
  an entropy of 1 corresponding to a disordered and random
  data set with all {eigenfeature, eigenassay}-pairs
  equally expressed} \item{\code{apply}:}{\code{character}
  containing whether the eigensystem should be computed for
  the actual data or the variance in the data}
  \item{\code{excludeEigenfeatures}:}{\code{numeric} vector
  containing eigenfeature 1 and 2 in case they capture
  >85\% of the data with eigenfeature 2 capturing at least
  15\%, otherwise numeric value containing eigenfeature 1}
  \item{\code{colorIdFeatures}:}{\code{numeric} vector or
  \code{factor} containing annotation information on the
  features} }
}

\section{Accessors}{
  \describe{ \item{}{\code{matrix(x)}, \code{matrix(x) <-
  value}} \item{}{\code{signMatrix(x)}, \code{signMatrix(x)
  <- value}} \item{}{\code{assayMatrix(x)},
  \code{assayMatrix(x) <- value}}
  \item{}{\code{featureMatrix(x)}, \code{featureMatrix(x)
  <- value}} \item{}{\code{eigenassays(x)},
  \code{eigenassays(x) <- value}}
  \item{}{\code{eigenexpressions(x)},
  \code{eigenexpressions(x) <- value}}
  \item{}{\code{eigenfeatures(x)}, \code{eigenfeatures(x)
  <- value}} \item{}{\code{assaycorrelations(x)},
  \code{assaycorrelations(x) <- value}}
  \item{}{\code{featurecorrelations(x)},
  \code{featurecorrelations(x) <- value}}
  \item{}{\code{fractions(x)}, \code{fractions(x) <-
  value}} \item{}{\code{entropy(x)}, \code{entropy(x) <-
  value}} \item{}{\code{apply(x)}, \code{apply(x) <-
  value}} \item{}{\code{excludeEigenfeatures(x)},
  \code{excludeEigenfeatures(x) <- value}}
  \item{}{\code{colorIdFeatures(x)},
  \code{colorIdFeatures(x) <- value}} }
}
\examples{
## Metabolomics starvation data obtained from http://genomics-pubs.princeton.edu/StarvationMetabolomics/Download.shtml
data(StarvationData)

## An object from class Eigensystem is obtained with the compute method
eigensystem <- compute(StarvationData)

## Obtain entropy
entropy(eigensystem)
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
\link{compute,Eigensystem-method}
}
\keyword{classes}
\keyword{methods}

