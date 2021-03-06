% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create-network.R
\name{transformNetwork}
\alias{transformNetwork}
\title{Add TPBN parameters to a BoolNet network.}
\usage{
transformNetwork(net, p = 0.01)
}
\arguments{
\item{net}{A network generated by \code{\link{createNetwork}}.}

\item{p}{Noise strength.}
}
\value{
A threshold PBN.
}
\description{
Add perturbation probability, threshold parameters, and Boolean functions
to the network, to turn it into a TPBN.
}
\examples{
# Create a BoolNet network and transform it.
net <- generateRandomNKNetwork(5, 2, "homogeneous")
net <- transformNetwork(net)
}
\seealso{
\code{\link{createNetwork}} \code{\link{A}} \code{\link{aToFunction}}
}
