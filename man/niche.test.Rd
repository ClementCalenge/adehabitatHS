\name{niche.test}
\alias{niche.test}
\title{Monte-Carlo Test on Parameters of the Ecological Niche}
\description{
  \code{niche.test} tests for the significance of two parameters of the
  ecological niche of a species (marginality and tolerance), using
  Monte-Carlo methods.  This is a bivariate test.
}
\usage{
niche.test(x, pts, nrep = 999, o.include = TRUE, \dots)
}
\arguments{
  \item{x}{a raster map of class \code{SpatialPixelsDataFrame}}
  \item{pts}{an object inheriting the class \code{SpatialPoints}}
  \item{nrep}{the number of permutations}
  \item{o.include}{logical, passed to \code{biv.test}.  If \code{TRUE},
    the origin is included in the plot}
  \item{\dots}{further arguments passed to \code{biv.test}}
}
\details{
  \code{niche.test} tests the significance of two parameters describing
  the ecological niche: the marginality
  (squared length of the vector linking the average available habitat
  conditions to the average used habitat conditions in the ecological
  space defined by the habitat variables), and the tolerance (inertia of
  the niche in the ecological space, i.e. the sum over all variables of
  the variance of used pixels). 
  
  At each step of the randomisation procedure, 
  the test randomly allocates the n points in the pixels of the map.
  The marginality and the tolerance are then recomputed on this
  randomised data set.
  
  Actual values are compared to random values with the help of the
  function \code{biv.test}.
}
\value{
Returns a list containing the following
components: 
\item{dfxy}{a data frame with the randomized values of marginality
  (first column) and tolerance (second column).}
\item{obs}{the actual value of marginality and tolerance.}
}
\author{
  Mathieu Basille \email{basille@ase-research.org}\cr
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\section{Warning}{
 \code{biv.test} uses the function \code{kde2d} of the package \code{MASS}.
}
\seealso{
  \code{\link{biv.test}} for more details on bivariate tests. 
  \code{\link{histniche}} for the histograms of the variables of the niche.
}
\examples{
\dontrun{
data(chamois)

niche <- niche.test(chamois$map,
                    chamois$locs,
                    side = "bottom")
names(niche)
}
}
\keyword{multivariate}
