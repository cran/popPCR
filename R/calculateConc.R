#' Target copies estimation
#'
#' Mean target copies per partition (lambda) is derived using Poisson distribution as lambda = -ln(nneg / ntot). Target copies in sample is then calculated as conc = lambda * volSamp/(volDrp * 1000).
#'
#' @param nneg numeric, negative droplet count
#' @param ntotal numeric, total droplet count
#' @param volSamp numeric, sample volume in microliter
#' @param volDrp  numeric, droplet (or partition) volume in nanoliter
#' @export
#' @return Returns a list with 2 named items `lambda` and `conc`
#' \itemize{
#'   \item lambda - numeric, vector of mean target copies per partition (lambda) and its lower and upper 95% confidence interval
#'   \item conc - numeric, vector of target copies in sample (based on the given sample volume (`volSamp`) and droplet volume (`volDrp`)) and its lower and upper 95% confidence interval
#' }
#' @examples
#' estimates <- calculateConc(5000, 20000, volSamp = 20, volDrp = 0.85)
#' estimates
#' #    Output:
#' #       $lambda
#' #          lambda     lower      upper
#' #       1.386294   1.362289   1.410299
#' #
#' #       $conc
#' #           conc      lower      upper
#' #       32618.69   32053.87   33183.51
calculateConc <- function(nneg, ntotal, volSamp, volDrp){
  res <- list()

  lambda <- -log(nneg/ntotal)
  lambda_lo <- lambda - 1.96 * sqrt((ntotal - nneg)/(ntotal * nneg))
  lambda_up <- lambda + 1.96 * sqrt((ntotal - nneg)/(ntotal * nneg))
  res$lambda <- c(lambda, lambda_lo, lambda_up)
  names(res$lambda) <- c("lambda", "lower", "upper")

  conc <- lambda * volSamp/volDrp * 1000
  conc_lo <- lambda_lo * volSamp/volDrp * 1000
  conc_up <- lambda_up * volSamp/volDrp * 1000
  res$conc <- c(conc, conc_lo, conc_up)
  names(res$conc) <- c("conc", "lower", "upper")
  return(res)
}
