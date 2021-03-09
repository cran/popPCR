#' EM Mixture Model fitting of dPCR droplet fluorescence
#'
#' Estimates target concentration by counting positive droplets with Poisson correction. Positive, negative, and rain populations are fitted using EM. Droplets are then classified using Maximum A Posteriori rule
#' @param x numeric, vector of droplet fluorescence amplitude
#' @param dist character, distribution of the mixture models ("normal", "skewed-normal", "t", "skewed-t")
#' @param volSamp numeric, sample volume in microliter
#' @param volDrp  numeric, droplet (or partition) volume in nanoliter
#' @param maxComponents numeric, maximum number of components (e.g. populations)
#' @param negProbThres numeric, if only one population was detected, then its assumed as a negative population. Droplets will be classified as positive if its probability given the population < negProbThres.
#' @param useOnlyNegProbThres logical, if TRUE, then droplets will be classified as positive if its probability given the leftmost population < negProbThres. Default is FALSE, i.e. classification is done by Maximum A Posteriori rule.
#' @return Returns a `result.popPCR` S4 class object with attributes
#' \itemize{
#'   \item classification - character, vector of droplet classification
#'   \item dist - character, user-specified parameter for the mixture model
#'   \item dropletCount - list, droplet classification count
#'   \item em - list, returned value of EMMIXskew's EmSkew()
#'   \item estConc - list, estimated target concentration as lambda and sample concentration (with 95% CI)
#'   \item G - numeric, number of components fitted
#'   \item memberProb - list, component membership probability of all droplets
#' }
#' @importFrom graphics plot
#' @importFrom methods new
#' @importFrom mvtnorm rmvnorm
#' @importFrom grDevices colorRampPalette contourLines densCols dev.new dev.off png
#' @importFrom graphics contour lines pairs par points smoothScatter
#' @importFrom stats as.dist cor cutree density hclust kmeans quantile rgamma rnorm uniroot
#' @useDynLib popPCR, .registration = TRUE
#' @export
#' @examples
#' library(popPCR)
#'
#' # Plot histograms of available data
#' hist(x_onePop, breaks = 100)
#' hist(x_twoPop, breaks = 100)
#' hist(x_multiPop, breaks = 100)
#'
#' # ---- Mixture model fitting ---- #
#' # Example 1. One population sample
#' result <- popPCR(x_onePop, dist = "t")
#' printSummaryConc(result)
#' #    Output:
#' #        Populations detected : 1
#' #        Total droplets : 8000
#' #        Positive : 1 (0.01%)
#' #        Negative : 7999 (99.99%)
#' #
#' #        Target copies in sample          : 2.9414 ( 95% CI: [ -2.8237 , 8.7064 ] )
#' #        Mean target copies per partition : 1e-04 ( 95% CI: [ -1e-04 , 4e-04 ] )
#' printSummaryFit(result)
#' #    Output:
#' #        Results of fitting a 1-component t mixture model
#' #
#' #        Negative Population
#' #        Mix prop. : 1
#' #        Mu        : 1024.1614
#' #        Sigma     : 35253.1747
#' #        Dof       : 2.005
#'
#' # (Option) increase negProbThres to classify negative droplets more strictly
#' result <- popPCR(x_onePop, dist = "t", negProbThres = 1e-4)
#' printSummaryConc(result)
#' #    Output:
#' #        Populations detected : 1
#' #        Total droplets : 8000
#' #        Positive : 691 (8.64%)
#' #        Negative : 7309 (91.36%)
#' #
#' #        Target copies in sample          : 2125.5312 ( 95% CI: [ 1966.9936 , 2284.0688 ] )
#' #        Mean target copies per partition : 0.0903 ( 95% CI: [ 0.0836 , 0.0971 ] )
#'
#' # Example 2. Two population sample
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
#'
#'
#' # Example 3. Multiple population sample
#' result <- popPCR(x_multiPop, dist = "t", maxComponents = 4)
#' printSummaryConc(result)
#' #    Output:
#' #        Populations detected : 4
#' #        Total droplets : 1814
#' #        Positive : 44 (2.43%)
#' #        Negative : 1252 (69.02%)
#' #        Rain (1) : 258 (14.22%)
#' #        Rain (2) : 260 (14.33%)
#' #
#' #        Target copies in sample          : 8724.5195 ( 95% CI: [ 7999.0578 , 9449.9812 ] )
#' #        Mean target copies per partition : 0.3708 ( 95% CI: [ 0.34 , 0.4016 ] )
#'
#' # In the output above, we see 2 rain populations! Let's examine its plot.
#' plot(stats::density(x_multiPop))
#' # We can see that Rain (1) is very close to the Negative population.
#' # Let's include droplets in Rain (1) in the negative droplet count.
#' nNegative <- result@dropletCount$neg + result@dropletCount$rain1
#' nTotal <- result@dropletCount$total
#' # Re-estimate concentration as follows
#' newEstimates <- calculateConc(nNegative, nTotal, volSamp = 20, volDrp = 0.85)
#' newEstimates
#' #    Output:
#' #       $lambda
#' #          lambda     lower     upper
#' #       0.1834247 0.1627763 0.2040731
#' #
#' #       $conc
#' #           conc    lower    upper
#' #       4315.875 3830.031 4801.719
#'
popPCR <- function(x, dist, volSamp=20, volDrp=0.85, maxComponents=Inf, negProbThres=0.0000001, useOnlyNegProbThres=FALSE){
  if(!dist %in% names(..distr)){
    stop(paste0("Invalid dist. Select one `", paste0(names(..distr), collapse="`, `"), "`"))
  }
  x <- stats::na.omit(x)
  x <- sort(x)

  modes <- ..getInitMus(x, maxComponents = maxComponents)
  G <- length(modes)
  initEMSkew <- list(mu = modes,
                     sigma = array(rep(1000, G), c(1,1,G)),
                     pro = rep(1/G, G),
                     dof = rep(30, G),
                     delta = t(matrix(rep(0,G))))
  df_x <- data.frame(Fluorescence=x)
  distr <- ..distr[[dist]]
  em <- EmSkew(df_x, init = initEMSkew, g = G, distr = distr, ncov = 3, debug = FALSE)

  probs <- ..getClusMemberProb(x, em, useOnlyNegProbThres, negProbThres)
  negMemberProb   <- probs$negMemberProb
  posMemberProb   <- probs$posMemberProb
  rainMemberProb  <- probs$rainMemberProb
  classification  <- probs$classification

  nneg <- sum(classification == "neg")
  npos <- sum(classification == "pos")
  dropletCount <- list(total = length(x), pos = npos, neg = nneg)
  for(rain in levels(classification)[grepl("rain", levels(classification))]){
    dropletCount[[rain]] <- sum(classification == rain)
  }

  em$mu    <- apply(em$mu,2,function(x){x})
  em$sigma <- apply(em$sigma,3,function(x){x})
  if("delta" %in% names(em)) em$delta <- apply(em$delta,2,function(x){x})

  r <- result.popPCR(estConc = calculateConc(nneg, length(x), volSamp, volDrp),
                     memberProb = list(droplet = x, negProb = negMemberProb, posProb = posMemberProb, rainProb = rainMemberProb),
                     G = G,
                     dist = dist,
                     em = em,
                     dropletCount = dropletCount,
                     classification = classification)

  printSummaryConc(r)

  invisible(r)
}

..getInitMus <- function(x, maxComponents=Inf){
  kerbw <- stats::bw.nrd0(x)
  if(kerbw < 50){
    kerbw <- 50
  }
  krn <- stats::density(x, bw = kerbw)
  krn <- rbind(krn$y, krn$x)

  piek <- ..findpeaks(krn[1, ], bw=20)$max.X
  piek <- rbind(piek, krn[1, piek])    #add peak heights
  piek <- rbind(piek, krn[2, piek[1,]])#add peak x-locations
  piek <- piek[, !(piek[2,] < max(piek[2,])/100)]
  if(is.null(dim(piek))){
    piek <- as.matrix(piek)
  }
  if(maxComponents < Inf){
    keep <- order(piek[2,], decreasing = TRUE)[1:maxComponents]
    keep <- keep[!is.na(keep)]
    piek <- piek[3,keep]
  }else{
    piek <- piek[3,]
  }
  return(piek)
}

..findpeaks <- function(vec, bw = 1, x.coo = c(1:length(vec))){
  # sourced from https://github.com/Gromgorgel/ddPCR/blob/92919f5c29f532d11d1f7074e621adf89350e839/Cloudy-V2-04.R

  #where bw = is box width, setting the sensitivity of the search
  ###set all vectors to null
  pos.x.max <- NULL ;	pos.y.max <- NULL ;	pos.x.min <- NULL ;	pos.y.min <- NULL
  ###Start of for loop:    we walk down the vector with a window of size "bw"
  for(i in 1:(length(vec)-1)){
    #check if we have reached the end of the vector
    if((i+1+bw)>length(vec)){sup.stop <- length(vec)}else{sup.stop <- i+1+bw}
    #check if we are at beginning of the vector
    if((i-bw) < 1){inf.stop <- 1}else{inf.stop <- i-bw}
    #select window in two parts: values beyond i (superior), and values before i (inferior)
    subset.sup <- vec[(i+1):sup.stop]
    subset.inf <- vec[inf.stop:(i-1)]
    ##############################################################
    #are ALL trailing data smaller than i?
    is.max   <- sum(subset.inf > vec[i]) == 0
    #are ALL leading data smaller than i?
    is.nomin <- sum(subset.sup > vec[i]) == 0
    #are ALL trailing data larger than i?
    no.max   <- sum(subset.inf > vec[i]) == length(subset.inf)
    #are ALL leading data larger than i?
    no.nomin <- sum(subset.sup > vec[i]) == length(subset.sup)
    ##############################################################
    #a maximum is found if  all data before and after i are smaller than i
    if(is.max & is.nomin){
      pos.x.max <- c(pos.x.max, x.coo[i])
      pos.y.max <- c(pos.y.max, vec[i])
    }
    #a maximum is found if  all data before and after i are larger than i
    if(no.max & no.nomin){
      pos.x.min <- c(pos.x.min, x.coo[i])
      pos.y.min <- c(pos.y.min, vec[i])
    }
  }#end of for loop
  ###Output
  return(list("max.X" = pos.x.max, "max.Y" = pos.y.max, "min.X" = pos.x.min, "min.Y" = pos.y.min))
}

..getClusMemberProb <- function(x, em, useOnlyNegProbThres, negProbThres){
  if(length(em$mu) > 1 && !useOnlyNegProbThres){
    clusterMem <- em$tau
    means <- em$modpts
    posClust <- which.max(means)
    negClust <- which.min(means)
    rainClust <- sort((1:length(means))[-c(posClust, negClust)])

    clusts <- c(posClust, negClust, rainClust)
    if(length(rainClust) == 0){
      names(clusts) <- c("pos", "neg")
    }else{
      names(clusts) <- c("pos", "neg", paste0("rain", 1:length(rainClust)))
    }
    clusts <- sort(clusts)
    classification <- factor(names(clusts[em$clust]), levels = names(clusts))

    # Correction for possible negative members beyond positive population
    a <- ifelse(classification=="neg",-1,1)
    b <- a[1:(length(a)-1)] + a[2:length(a)]
    erroneousClust <- which(b == 0)[2]
    if(!is.na(erroneousClust)){
      temp <- classification[erroneousClust:length(classification)]
      temp[temp=="neg"] <- "pos"
      classification[erroneousClust:length(classification)] <- temp
    }

    posMemberProb <- clusterMem[,posClust]
    negMemberProb <- clusterMem[,negClust]
    rainMemberProb <- clusterMem[,rainClust]
  }else{
    # If population is 1, assumes that the population detected is the negative distribution
    mean <- em$mu[1,1]
    sigma <- em$sigma[1,1,1]
    df <- em$dof[1]
    x_mayBePos <- x[x > mean]
    if(em$distr == "mvn"){
      d <- ddmvn(matrix(x_mayBePos, ncol=1), n=length(x_mayBePos), p=1, mean=mean, cov=sigma)
    }
    else if(em$distr == "msn"){
      d <- ddmsn(matrix(x_mayBePos, ncol=1), n=length(x_mayBePos), p=1, mean=mean, cov=sigma, del =  em$delta[1,1])
    }
    else if(em$distr == "mvt"){
      d <- ddmvt(matrix(x_mayBePos, ncol=1), n=length(x_mayBePos), p=1, mean=mean, cov=sigma, nu = df)
    }
    else if(em$distr == "mst"){
      d <- ddmst(matrix(x_mayBePos, ncol=1), n=length(x_mayBePos), p=1, mean=mean, cov=sigma, nu = df, del =  em$delta[1,1])
    }
    surenegs <- length(x)-length(x_mayBePos)
    i_pos <- which(d < negProbThres) + surenegs
    i_neg <- setdiff(1:length(x), i_pos)

    negMemberProb <- c(rep(1, surenegs), d)
    posMemberProb <- 1 - negMemberProb
    rainMemberProb <- NA

    classification <- character(length(x))
    classification[i_neg] <- "neg"
    classification[i_pos] <- "pos"
    classification <- factor(classification, levels = c("pos", "neg"))
  }
  return(list(negMemberProb = negMemberProb, posMemberProb = posMemberProb, rainMemberProb = rainMemberProb, classification = classification))
}

..distr <- c("normal"        = "mvn",
             "skewed-normal" = "msn",
             "t"             = "mvt",
             "skewed-t"      = "mst")

result.popPCR <- methods::setClass("result.popPCR", slots=list(
  classification = "factor",
  dist = "character",
  dropletCount = "list",
  em = "list",
  estConc = "list",
  G = "numeric",
  memberProb = "list"
))
