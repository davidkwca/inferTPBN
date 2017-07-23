##' Infer a Threshold PBN from a list of timeseries.
##'
##'
##' This is a wrapper for the whole inference procedure and makes use of
##' \code{\link{innovateGeneUntilSaturated}}.
##' For every gene it infers regulatory functions until the upper bound is
##' reached or there is no more improvement from adding more functions.
##'
##' @param ts.multi A list of timeseries.
##' @param p The noise probability to use for the network.
##' @param L.max The maximum number of regulatory networks to try to infer.
##' @param n.cores The number of cores to be used.
##' @param seed The random seed to use.
##' @param partial Whether or not to *also* infer a network where only partial
##' optimization over the last threshold parameter is performed. Defaults to
##'   FALSE.
##' @param verbal Show progress report? Defaults to FALSE.
##' @return A list containing the inferred network and the time taken to perform
##'   the inference.
##' @seealso \code{\link{FullRun}}
##' @examples
##' net.true <- createNetwork(inputProbabilities=c(0.5, 0.5), n=5, k=2)
##' ts.multi <- simulateNetwork(net.true, 10, 20)
##' inferred.list <- inferPBN(ts.multi, n.cores=1)
##' net.inferred <- inferred.list$net.inferred
##' time <- inferred.list$time.complete
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
  ## print(net.inferred$interactions[[1]])

  registerDoParallel(n.cores)
  time.complete <-
    system.time({
      temp <- foreach(g = 1:n) %dopar% {
        innovateGeneUntilSaturated(ts.multi, net.inferred, g,
                                   partial=FALSE, verbal=verbal)
      }
      improvement <- rep(0, n)
      for (g in 1:n){
        if (temp[[g]]$interactions[[g]][[1]]$input != 0){
          improvement[g] <- 1
        }
        net.inferred$interactions[[g]] <- temp[[g]]$interactions[[g]]
      }
      ## Only up to L.max new regulatory functions are inferred per gene.
      for (L in 1:L.max){
        ## New regulatory functions are inferred only where improvements have
        ## been made on the previous round
        improved <- (1:n)[improvement > 0]
        temp <- foreach(g = improved) %dopar% {
          innovateGeneFunction(ts.multi, net.inferred, g,
                               partial=FALSE, verbal=verbal)
        }
        improvement.temp <- rep(0, n)
        for (j in 1:length(improved)){
          g <- improved[j]
          if (length(net.inferred$interactions[[g]])
              < length(temp[[j]]$interactions[[g]])){
            improvement.temp[g] <- 1
          }
          net.inferred$interactions[[g]] <- temp[[j]]$interactions[[g]]
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

  return(returns)
}
