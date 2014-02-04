\name{StarvationData}
\docType{data}
\alias{StarvationData}
\title{Human HeLa Cell Cycle Expression}
\usage{
  data(StarvationData)
  StarvationData
}
\description{
Metabolic response to starvation was studied in two microbes, Escherichia coli and Saccharomyces cerevisae, to determine whether metabolome response to nutrient deprivation is similar across both organisms. Sixty-eight cellular metabolites were analyzed by LC-MS/MS in both bacteria and yeast, after nutrient starvation with carbon and nitrogen. Cells were sampled for 8 hours.
}
\format{
	An ExpressionSet with metabolite x sample intensity data and sample information. Intensity data comprise the logtransformed relative metabolite concentration changes with respect to experiment initiation at time point 0 hrs. For the samples, information is provided on species (Bacteria vs Yeast), starvation (Carbon vs Nitrogen) and time point in hrs.
}
\source{
	Brauer MJ, Yuan J, Bennett BD, et al. Conservation of the metabolomic response to starvation across two divergent microbes. Proc Nat Acad Sci U.S.A. 103(51), 19302-19307 (2006).
}
\keyword{datasets}
\author{
  Anneleen Daemen, Matthew Brauer
}
