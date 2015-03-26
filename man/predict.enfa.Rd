\name{predict.enfa}
\alias{predict.enfa}
\title{Habitat Suitability Maps Built from the ENFA}
\description{
  \code{predict.enfa} computes habitat suitability maps using the
  Ecological-Niche Factor Analysis and the Mahalanobis distances
  method.
}
\usage{
\method{predict}{enfa}(object, map, nf, \dots)
}
\arguments{
  \item{object}{an object of class \code{enfa}}
  \item{map}{an object of class \code{SpatialPixelsDataFrame}.} 
  \item{nf}{the number of axes of specialization kept for the
    predictions.  By default, all axes kept in \code{object} are
    used}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{

  The predictions are based on the position of the niche defined by the
  ENFA within the multidimensional space of environmental variables.  The
  ENFA produces row coordinates for each pixel, which are used with the
  function \code{mahalanobis}.  For each pixel, this function computes the
  Mahalanobis distances from the barycentre of the niche.
  
  Actually, the function \code{predict.enfa} is identical to the function
  \code{mahasuhab}, except that the habitat suitability map is computed
  using the axes of the ENFA, instead of the raw data.

  Note that the MADIFA allows a more consistent factorial decomposition of
  the Mahalanobis distances.
}
\value{
Returns a raster map of class \code{SpatialPixelsDataFrame}.
}
\references{
  Clark, J.D., Dunn, J.E. and Smith, K.G. (1993) A multivariate model of
  female black bear habitat use for a geographic information
  system. \emph{Journal of Wildlife Management}, \bold{57}, 519--526.

  Hirzel, A.H., Hausser, J., Chessel, D. & Perrin, N. (2002)
  Ecological-niche factor analysis: How to compute habitat-suitability
  maps without absence data? \emph{Ecology}, \bold{83}, 2027--2036.  
}
\author{Mathieu Basille \email{basille@ase-research.org}}
\seealso{
\code{\link[stats]{mahalanobis}} for information on the computation of
Mahalanobis distances.  \code{\link{mahasuhab}} for more details on the
computation of habitat suitability maps using the Mahalanobis distances.
\code{\link{madifa}} for a more consistent factorial decomposition of
the Mahalanobis distances
}
\examples{
\dontrun{
data(lynxjura)

map <- lynxjura$map

## We keep only "wild" indices.
locs <- lynxjura$locs
locs <- locs[slot(locs, "data")[,2]!="D",]

pr <- slot(count.points(locs, map), "data")[,1]

(enfa1 <- enfa(dudi.pca(slot(map, "data"), scannf=FALSE),
               pr, scannf = FALSE))

## Compute the prediction
pred <- predict(enfa1, map)
image(pred)
contour(pred, col="green", add=TRUE)
points(locs, pch = 3)

}
}
\keyword{multivariate}
