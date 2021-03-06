% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateConc.R
\name{calculateConc}
\alias{calculateConc}
\title{Target copies estimation}
\usage{
calculateConc(nneg, ntotal, volSamp, volDrp)
}
\arguments{
\item{nneg}{numeric, negative droplet count}

\item{ntotal}{numeric, total droplet count}

\item{volSamp}{numeric, sample volume in microliter}

\item{volDrp}{numeric, droplet (or partition) volume in nanoliter}
}
\value{
Returns a list with 2 named items \code{lambda} and \code{conc}
\itemize{
\item lambda - numeric, vector of mean target copies per partition (lambda) and its lower and upper 95\% confidence interval
\item conc - numeric, vector of target copies in sample (based on the given sample volume (\code{volSamp}) and droplet volume (\code{volDrp})) and its lower and upper 95\% confidence interval
}
}
\description{
Mean target copies per partition (lambda) is derived using Poisson distribution as lambda = -ln(nneg / ntot). Target copies in sample is then calculated as conc = lambda * volSamp/(volDrp * 1000).
}
\examples{
estimates <- calculateConc(5000, 20000, volSamp = 20, volDrp = 0.85)
estimates
#    Output:
#       $lambda
#          lambda     lower      upper
#       1.386294   1.362289   1.410299
#
#       $conc
#           conc      lower      upper
#       32618.69   32053.87   33183.51
}
