\name{domain}
\alias{domain}
\title{Estimation of the Potential Distribution of a Species }
\description{
\code{domain} uses the DOMAIN algorithm to estimate the potential
distribution of a species based on a list of species occurrences and on
maps of the area.
}
\usage{
domain(x, pts, type = c("value", "potential"), thresh = 0.95)
}
\arguments{
  \item{x}{an object of class \code{SpatialPixelsDataFrame}}
  \item{pts}{a data frame giving the x and y coordinates of the species
    occurrences.}
  \item{type}{a character string. The \code{"value"} of the suitability
    may be returned 
    or the \code{"potential"} area of distribution}
  \item{thresh}{if \code{value = "potential"}, a threshold value should be
    supplied for the suitability (by default 0.95)}
}
\details{
  This function implements the DOMAIN algorithm described in Carpenter
  et al. (1993).
}
\value{
  Returns an object of class \code{SpatialPixelsDataFrame}.
}
\section{Warning }{
  \code{domain} is restricted to maps
  containing only numerical variables (i.e. no factors).
}
\references{
  Carpenter, G., Gillison, A.N. and Winter, J. (1993)
  DOMAIN: a flexible modelling
  procedure for mapping potential distributions of plants and
  animals. \emph{Biodiversity and conservation}, \bold{2}, 667--680.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link{mahasuhab}}   }

\examples{

## Preparation of the data
data(lynxjura)
map <- lynxjura$map
pts <- lynxjura$locs

## View of the data
image(map)
title(main="Elevation")
points(pts, pch = 3)

## Estimation of habitat suitability map
hsm <- domain(map, pts)

image(hsm, col = grey((1:256)/256))
contour(hsm, add = TRUE)

## Lighter areas are the most used areas

## Potential distribution
hsm <- domain(map, pts, type = "potential", thresh = 0.98)
image(hsm, col = "orange")
title(main = "Habitat suitability map")
points(pts, pch = 3)

}
\keyword{spatial}
\keyword{multivariate}
