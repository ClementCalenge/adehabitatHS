\name{mahasuhab}
\alias{mahasuhab}
\title{Habitat Suitability Mapping with Mahalanobis Distances.}
\description{
  This function computes the habitat suitability map of an area
  for a species, given a set of locations of the species occurences (Clark
  et al. 1993).  This function is
  to be used in habitat selection studies, when animals are not
  identified.
}
\usage{
mahasuhab(x, pts, type = c("distance", "probability"))
}
\arguments{
  \item{x}{a raster map of class \code{SpatialPixelsDataFrame}}
  \item{pts}{a data frame with two columns, giving the coordinates
    of the species locations}
  \item{type}{a character string. Whether the raw \code{"distance"} should be
    returned, or rather the \code{"probability"} (see details). }
}
\details{
  Let assume that a set of locations of the species on an area is available
  (gathered on transects, or during the monitoring of the population,
  etc.).  If we assume that the probability of detecting an individual
  is independent from the habitat variables, then we can consider that
  the habitat found at these sites reflects the habitat use by the animals.
  
  The Mahalanobis distances method has become more and more popular
  during the past few years to derive habitat suitability maps.  The
  niche of a species is defined as the probability density function of
  presence of a species in the multidimensionnal space defined by the
  habitat variables.  If this function can be assumed to
  be multivariate normal, then the mean vector of this distribution
  corresponds to the optimum for the species.

  The function \code{mahasuhab} first computes this mean vector as well
  as the variance-covariance matrix of the niche density function, based
  on the value of habitat variables in the sample of locations.
  Then, the *squared* Mahalanobis distance from this optimum is computed
  for each pixel of the map.  Thus, the smaller this squared
  distance is for a given pixel, and the better is the habitat in this
  pixel.

  Assuming multivariate normality, squared Mahalanobis distances are
  approximately distributed as Chi-square with n degrees of freedom,
  where n equals the number of habitat characteristics (see the section
  note below on this question).  If the
  argument \code{type = "probability"}, maps of these p-values are
  returned by the function. As such these are the probabilities of a
  larger squared Mahalanobis distance than that observed when x is
  sampled from the niche.
}
\value{
  Returns a raster map of class \code{SpatialPixelsDataFrame}.
}
\references{
  Clark, J.D., Dunn, J.E. and Smith, K.G. (1993) A multivariate model of
  female black bear habitat use for a geographic information
  system. \emph{Journal of Wildlife Management}, \bold{57}, 519--526.

  Etherington, T.R. (2019) Mahalanobis distances and ecological niche
  modelling: correcting a chi-squared probability error. \emph{Peer J},
  \bold{7}, e6678.
}
\author{ Clement Calenge \email{clement.calenge@ofb.gouv.fr} }
\note{

  The computation of the squared Mahalanobis distances inverts the
  variance-covariance matrix of the niche density function (see
  \code{?mahalanobis}).
  It is therefore important that the habitat variables considered
  are not too correlated among each other.  When the habitat variables
  are too correlated, the variance-covariance matrix is singular and
  cannot be inverted.

  Note also that it is recommended to scale the variables before the
  computation, so that they all have the same variance, and therefore
  the same weight in the analysis (see examples below).
 
  Finally, note that in versions of adehabitatHS prior to 0.3.16,
  \code{mahasuhab} incorrectly calculated the probability, when
  \code{type = "probability"}, based on a Chi-square with (n-1) degrees
  of freedom instead of n degrees of freedom (see Etherington 2019 on
  this issue). 
}
\seealso{ \code{\link{madifa}} and \code{\link{dunnfa}} for factor
  analyses of the Mahalanobis distances, \code{\link{domain}} for
  another method of habitat suitability mapping,
  \code{\link[stats]{mahalanobis}} for
  information on the computation of Mahalanobis distances. }
\examples{
## loads the data
data(lynxjura)
ma <- lynxjura$map
lo <- lynxjura$locs[,1:2]

## We first scale the maps
slot(ma, "data") <- dudi.pca(slot(ma, "data"), scannf=FALSE)$tab

## habitat suitability mapping
hsm <- mahasuhab(ma, lo, type = "probability")
image(hsm)
title(main = "Habitat suitability map for the Lynx")
points(lo, pch = 3)

}
\keyword{spatial}
\keyword{multivariate}
