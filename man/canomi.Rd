\name{canomi}
\alias{canomi}
\alias{print.canomi}
\alias{plot.canomi}
\title{
  Canonical OMI analysis
}
\description{
  This function performs a canonical OMI analysis (outlying mean index).
}
\usage{
canomi(dudiX, Y, scannf = TRUE, nf = 2)
\method{print}{canomi}(x, \dots)
\method{plot}{canomi} (x, xax = 1, yax = 2, \dots)
}
\arguments{
  \item{dudiX}{an object of class \code{dudi}}
  \item{Y}{a a data frame Resource units-animals according to
    \code{dudiX$tab} with no columns of zero}
  \item{scannf}{a logical value indicating whether the eigenvalues bar plot
    should be displayed}
  \item{nf}{if scannf FALSE, an integer indicating the number of kept
    axes}
  \item{x}{an object of class \code{canomi} returned by the function
    \code{canomi}}
  \item{xax}{the number of the x-axis}
  \item{yax}{the number of the x-axis}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  The canonical OMI analysis is similar to the function \code{niche},
  from the package \code{ade4}.  The principle of this analysis is the
  following.  A set of N resource units (RUs) are available to the
  animals of the study.  Each resource unit is described by P
  environmental variables.  Therefore, the N resource units define a
  cloud of N points in the P-dimensionnal space defined by the P
  variables. We call this space "ecological space".

  Moreover, the use of the N resource units is known (or sampled) for
  a sample of K animals (e.g., using radio-tracking).  These utilization
  weights for each RU (rows) and each animal (column) define a table
  \code{Y}.  For a given animal, the set of resource units used define
  the "niche" of the animal.  The vector connecting the centroid (mean)
  of the available RUs to the centroid of the RUs used by this animal is
  named "marginality vector" (and its squared length is named
  "marginality" or "outlying mean index").

  The \code{canomi} first distorts the ecological space so that the
  available resource units take a standard spherical shape (by first
  performing a principal component analysis). Then, in this distorted
  space, a non-centred principal component analysis of the marginality
  vectors is performed. The canonical OMI analysis finds the directions
  in the distorted ecological space where the marginality is, in
  average, the largest.
}
\value{
  \code{canomi} returns a list of the class \code{canomi}, with the
  following components:
  \item{call}{original call.}
  \item{rank}{an integer indicating the rank of the studied matrix}
  \item{nf}{an integer indicating the number of kept axes}
  \item{eig}{a vector with all the eigenvalues of the analysis.}  
  \item{tab}{a data frame with n rows (n animals) and p columns (p
  environmental variables).}
  \item{li}{animals coordinates, data frame with n rows and nf columns.}
  \item{l1}{animals normed coordinates, data frame with n rows and nf columns.}
  \item{c1}{column scores, data frame with p rows and nf columns.}
  \item{cor}{the correlation between the canomi axes and the original
    variables}
  \item{ls}{a data frame with the resource units coordinates}
  \item{cm}{The variables metric used in the analysis (e.g.
    \code{ls = dudiX\$tab\%*\%cm\%*\%c1})}
  \item{as}{a data frame with the axis upon niche axis}
}
\references{
  Chessel, D. 2006. Calcul de l'outlier mean index. Consultation
  statistique avec le logiciel R.
}
\author{Clement Calenge \email{clement.calenge@ofb.gouv.fr} }
\seealso{
  \code{\link[ade4]{dudi}} for class \code{dudi},
  \code{\link[ade4]{niche}} for classical OMI analysis
}
\examples{
## The data
data(puech)
locs <- puech$relocations
maps <- puech$maps

## the maps
mimage(maps)

## the relocations of the wild boar:
image(maps)
points(locs, col=as.numeric(slot(locs, "data")[,1]))


## count the number of relocations
## in each pixel of the maps
cp <- count.points(locs, maps)

## gets the data:
dfavail <- slot(maps, "data")
dfused <- slot(cp, "data")

## a preliminary principal component analysis of the data:
dud <- dudi.pca(dfavail, scannf=FALSE)

## The analysis:
nic <- canomi(dud, dfused, scannf=FALSE)
nic

## Plot the results:
plot(nic)


}
\keyword{multivariate}
