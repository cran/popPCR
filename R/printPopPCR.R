#' Print result summary of popPCR
#'
#' Summarizes the number of populations detected, total droplets, and number of classified positive, negative, and rain droplets. Also calculates the target copies in sample and the mean target copies per partition (lambda).
#'
#' @param result.popPCR returned value of popPCR()
#' @export
#' @return NULL
#' @examples
#' result <- popPCR(x_twoPop, dist = "t")
#' printSummaryConc(result)
#' #    Output:
#' #        Populations detected : 2
#' #        Total droplets : 10254
#' #        Positive : 8693 (84.78%)
#' #        Negative : 1561 (15.22%)
#' #
#' #        Target copies in sample          : 44290.3819 ( 95% CI: [ 43215.6408 , 45365.1231 ] )
#' #        Mean target copies per partition : 1.8823 ( 95% CI: [ 1.8367 , 1.928 ] )
printSummaryConc <- function(result.popPCR){
  r <- result.popPCR
  posPct <- round(r@dropletCount$pos / r@dropletCount$total*100, 2)
  negPct <- round(r@dropletCount$neg / r@dropletCount$total*100, 2)

  lamb <- sapply(r@estConc$lambda, round, 4)
  conc <- sapply(r@estConc$conc, round, 4)

  summaryRain <- ""
  for(rain in levels(r@classification)[grepl("rain", levels(r@classification))]){
    i <- gsub("rain", "", rain)
    rainPct <- round(r@dropletCount[[rain]] / r@dropletCount$total*100, 2)
    summaryRain <- paste0(summaryRain, "Rain (", i, ") : ", sum(r@dropletCount[[rain]]), " (",rainPct,"%)\n      ")
  }

  cat(paste0("
      Populations detected : ", r@G, "
      Total droplets : ", r@dropletCount$total, "
      Positive : ", r@dropletCount$pos, " (",posPct,"%)
      Negative : ", r@dropletCount$neg, " (",negPct,"%)
      ",summaryRain,"
      Target copies in sample          : ", conc[['conc']],  " ( 95% CI: [ ", conc[['lower']]," , ", conc[['upper']]," ] )
      Mean target copies per partition : ", lamb[['lambda']]," ( 95% CI: [ ", lamb[['lower']]," , ", lamb[['upper']]," ] )"))
  writeLines("")
}

#' Print fitted mixture model estimates from popPCR
#'
#' Summarizes the number of populations fitted and their estimate distribution parameters. If only 1 population was detected, then it is assumed and is identified to be a negative population. If 2 populations were detected, then the leftmost is identified as the Negative Population and the rightmost is the Positive Population. If 3 or more populations were detected, then the populations between the leftmost and the rightmost will be considered as Rain Populations; which are numbered to make it identifiable in case of multiple Rain Populations (i.e. Rain (1) and Rain (2)).
#'
#' @param result.popPCR returned value of popPCR()
#' @export
#' @return NULL
#' @examples
#' result <- popPCR(x_twoPop, dist = "t")
#' printSummaryFit(result)
#' #    Output:
#' #        Results of fitting a 2-component t mixture model
#' #
#' #        Negative Population
#' #        Mix prop. : 0.1522
#' #        Mu        : 2136.7435
#' #        Sigma     : 4126.8357
#' #        Dof       : 12.3562
#' #
#' #        Positive Population
#' #        Mix prop. : 0.8478
#' #        Mu        : 7580.1275
#' #        Sigma     : 42621.1894
#' #        Dof       : 2.415
#' result <- popPCR(x_multiPop, dist = "t", maxComponents = 4)
#' printSummaryFit(result)
#' #     Output:
#' #        Results of fitting a 4-component t mixture model
#' #
#' #        Negative Population
#' #        Mix prop. : 0.6896
#' #        Mu        : 1452.1416
#' #        Sigma     : 12526.8931
#' #        Dof       : 21.3612
#' #
#' #        Rain (1) Population
#' #        Mix prop. : 0.142
#' #        Mu        : 2142.1118
#' #        Sigma     : 10762.5474
#' #        Dof       : 186.2947
#' #
#' #        Rain (2) Population
#' #        Mix prop. : 0.1457
#' #        Mu        : 5119.0039
#' #        Sigma     : 334959.2499
#' #        Dof       : 2.3626
#' #
#' #        Positive Population
#' #        Mix prop. : 0.0227
#' #        Mu        : 8505.9682
#' #        Sigma     : 192858.9044
#' #        Dof       : 149.8677
printSummaryFit <- function(result.popPCR){
  r <- result.popPCR

  cat(paste0("Results of fitting a ", r@G ,"-component ", r@dist, " mixture model"))

  hasDelta <- r@em$distr == "msn" || r@em$distr == "mst"
  hasDof   <- r@em$distr == "mvt" || r@em$distr == "mst"
  for(i in 1:r@G){
    if(i == 1){
      popn <- "Negative Population"
    }else if(i == r@G){
      popn <- "Positive Population"
    }else{
      popn <- paste0("Rain (", i-1, ") Population")
    }
    cat(paste0("

      ",popn,"
        Mix prop. : ", round(r@em$pro[i],4),"
        Mu        : ", round(r@em$mu[i],4),"
        Sigma     : ", round(r@em$sigma[i],4)))

    if(hasDof)   cat(paste0("\n        Dof       : ", round(r@em$dof[i],4)))
    if(hasDelta) cat(paste0("\n        Delta     : ", round(r@em$delta[i],4)))
  }
  writeLines("")
}
