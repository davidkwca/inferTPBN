##' Infer a Threshold PBN from a list of timeseries. Makes use of
##'   \code{\link{innovateGeneUntilSaturated}}
##'
##' det
##' @title inferPBN
##' @param ts.multi A list of timeseries.
##' @param p The noise probability to use for the network.
##' @param n.regulator.max
##' @param n.cores The number of cores to be used.
##' @param seed The random seed to use.
##' @param partial Whether or not to *also* infer a network where only partial
##' optimization over the last threshold parameter is performed. Defaults to
##'   FALSE.
##' @param verbal Show progress report? Defaults to FALSE.
##' @return A list containing the inferred network and the time taken to perform
##'   the inference.
##' @examples \dontrun{inferPBN(ts.multi)}
inferPBN <- function(ts.multi,
                     p=0.01,
                     L.max=2,
                     n.cores=detectCores()-1,
                     seed=111,
                     partial=FALSE, verbal=FALSE){
  returns <- list()
  set.seed(seed)
  n <- nrow(ts.multi[[1]])
  net.empty <- createNetwork(inputProbabilities=1, n=n , k=0,
                             topology="homogeneous",
                             noIrrelevantGenes = TRUE)

  for (g in 1:n){
    net.empty$interactions[[g]][[1]]$input <- 0
    net.empty$interactions[[g]][[1]]$expression <- NULL
  }

  net.inferred <- transformNetwork(net.empty)

  registerDoParallel(n.cores)
  time.complete <-
    system.time({
      temp <- foreach(g = 1:n) %dopar% {
        innovateGeneUntilSaturated(ts.multi, net.inferred, g,
                                   partial=FALSE, verbal=verbal)
      }
      improvement <- rep(0, n)
      for (g in 1:n){
        if (temp[[g]]$interactions[[g]]$input != 0){
          improvement[g] <- 1
        }
        net.inferred$interactions[[g]] <- temp[[g]]$interactions[[g]]
      }
      ## Only up to L.max new regulatory functions are inferred per gene.
      for (L in 1:L.max){
        ## New regulatory functions are inferred only where improvements have
        ## been made on the previous round
        temp <- foreach(g = (1:n)[improvement > 0]) %dopar% {
          innovateGeneFunction(ts.multi, net.inferred, g,
                               partial=FALSE, verbal=verbal)
        }
        improvement.temp <- rep(0, n)
        for (g in (1:n)[improvement > 0]){
          if (length(net.inferred$interactions[[g]])
              < length(temp[[g]]$interactions[[g]])):
               improvement.temp[g] <- 1
          net.inferred$interactions[[g]] <- temp[[g]]$interactions[[g]]
        }
        improvement <- improvement.temp
      }
    }
    )
  stopImplicitCluster()

  returns$net.inferred <- net.inferred
  returns$time.complete <- time.complete


  if(partial){
    net.inferred.partial <- transformNetwork(net.empty)

    time.partial <-
      system.time(
        for (g in 1:n){
          net.inferred.partial <-
            innovateGeneUntilSaturated(ts.multi,
                                       net.inferred.partial, g,
                                       partial=TRUE, verbal=verbal)
        })

    time.partial.update <-
      system.time(
        for (g in 1:n){
          net.inferred.partial <-
            innovateGeneFunction(ts.multi,
                                 net.inferred.partial, g,
                                 partial=TRUE, verbal=verbal)
        })
    returns$net.inferred.partial <- net.inferred.partial
    returns$time.complete.update.partial <- time.complete.update.partial
    returns$time.complete.partial <- time.complete.partial
  }

  ## returns$net.bestfit <- net.bestfit
  ## returns$time.bestfit <- time.bestfit

  return(returns)
}
