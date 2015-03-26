\name{squirreloc}
\alias{squirreloc}
\docType{data}
\title{Radio-tracking of squirrels}
\description{
  This data set contains the trajectories of 15 radio-monitored
  squirrels, as well as the vector maps of habitat composition.
}
\usage{data(squirreloc)}
\format{
  This data set is a list of two objects:
  \code{locs} is a \code{SpatialPointsDataFrame} containing the
  relocations of 15 squirrels
  
  \code{map} is an object of class \code{SpatialPolygonsDataFrame}
    containing the habitat composition of the area. The habitat types
    and colour coding are stored in the data.frame
}
\details{
  The dataset \code{squirreloc} comes from the Ranges VI software.  It
  has been used to illustrate the compositional analysis (see
  \code{?compana}) and the eigenanalysis of selection ratios (see
  \code{?eisera}).  See also the dataset \code{squirrel}.
}
\source{
  Kenward, R.E., South, A.B. and Walls, S.S. (2003). Ranges6 v1.2 : For
  the analysis of tracking and location data. Online manual. Anatrack
  Ltd. Wareham, UK. ISBN 0-9546327-0-2.
}
\keyword{datasets}
