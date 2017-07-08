##' Optimize a subset of threshold parameters.
##'
##'
##' While holding everything but the subset constant, optimize over this subset.
##' This essentially just tries every possible combination for the values on
##' this subset. However, zeroes are not allowed, as they have been found never
##' to come up while making the inference much, much, slower.
##'
##' @param ts.multi A list of timeseries.
##' @param net A network.
##' @param g A gene.
##' @param i The index of the specific regulatory function for the optimization.
##' @param subset The subset of parameters over which to optimize.
##' @return A vector containing the best threshold parameters found.
bestPartialA <- function(ts.multi, net, g, i, subset){
  gene <- net$genes[g]
  parents <- net$interactions[gene][[1]][[i]]$input
  parents <- parents[subset]
  c <- net$interactions[gene][[1]][[i]]$probability
  p.flip <- net$flipProbability
  p.perturbation <- 1 - (1 - p.flip)^length(net$genes)
  k <- length(parents)
  a.allowed <- c(-1, 1)

  ## I noticed that if the coefficient of a gene in regulating another one has
  ## been found not to be zero at an earlier iteration, later ones never changed
  ## this.
  ## However, the value tends to flip signs often enough that doing a "full"
  ## optimization is still a good thing.
  if(k == 1){
    a.allowed <- c(0, a.allowed)
  }

  A <- lapply(1:k, function(x) a.allowed)
  A <- expand.grid(A)

  ## Compute the loss function for each possible a and return minimum.
  ## If multiple produce the same error, the one with most zeros is returned due
  ## to the ordering used in generating A.
  m <- apply(A, 1,
             function(x){
               net$interactions[gene][[1]][[i]]$a[subset] <- x
               net$interactions[gene][[1]][[i]]$func <-
                 aToFunctionGene(net, gene, i)
               return(LossGeneTotal(ts.multi, net, g))
             }
             )
  a.min <- A[which.min(m),]
  return(as.numeric(a.min))
}


##' Jointly optimize all threshold parameters.
##'
##'
##' Wrapping \code{\link{bestPartialA}}, this doesn't hold any values constant.
##'
##' @param ts.multi A list of timeseries.
##' @param net A network.
##' @param g A gene.
##' @param i The index of the regulatory function.
##' @return A vector containing the best threshold parameters found.
bestA <- function(ts.multi, net, g, i){
  parents <- net$interactions[[g]][[i]]$input
  k <- length(parents)
  return(bestPartialA(ts.multi, net, g, i, 1:k))
}

##' Find optimal probabilities for a set of regulatory functions.
##'
##'
##' Using standard optimization \code{solnp}, this tries to find the
##' probabilities for a single gene which minimizes the loss function.
##'
##' @param ts.multi List of timeseries.
##' @param net A network.
##' @param g A gene.
##' @return A vector containing the probabilities optimizing the loss function.
bestC <- function(ts.multi, net, g){
  net.new <- net

  objective <- function(c){
    for (i in 1:length(c)){
      net.new$interactions[g][[1]][[i]]$probability <- c[i]
    }

    ## Entropy as regularizer because we would like to have as uneven
    ## probabilities as possible.
    loss <- LossGeneTotal(ts.multi, net.new, g)
  }

  ## The equality constraint that we want c to be a probability.
  equal <- function(c){
    sum(c)
  }

  n.regu <- length(net.new$interactions[[g]])
  c <- runif(n.regu, 0, 1)
  c <- c / sum(c)
  ## Don't show all the optimization outputs...
  capture.output(c.best
                 <- solnp(c, objective, equal, 1, LB=rep(0, n.regu), UB=rep(1, n.regu)))
  return(c.best)
}


##' Find best probabilities. Prune those with small values.
##'
##'
##' Using \code{\link{bestC}}, this finds the best probabilities. Afterwards, it
##' prunes all functions whose probabilities is below a user-defined threshold.
##'
##' @title Optimize probabilities for regulatory functions to be used.
##' @param ts.multi A list of timeseries.
##' @param net A network.
##' @param g A gene.
##' @param prune Whether or not to prune regulatory functions with low
##'   probabilities assigned to them.
##' @param thresh The threshold below which to prune functions if prune
##'   is TRUE.
##' @return A network whose probabilities have been optimized.
optimizeC <- function(ts.multi, net, g, prune=TRUE, thresh=0.1){
  ## If there is only a single regulatory set, no optimization needs to be done.
  if (length(net$interactions[g][[1]]) > 1){
    c.net <- bestC(ts.multi, net, g)$pars
    c.sum <- sum(c.net[c.net > thresh])
    for (i in 1:length(c.net)){
      net$interactions[g][[1]][[i]]$probability <- c.net[i]
    }

    ## Drop negligible functions.
    if(prune){
      q <- sapply(net$interactions[g][[1]], function(x) x$probability < thresh)
      net$interactions[g][[1]][q] <- NULL

      ## Renormalize probabilities.
      for (i in 1:length(net$interactions[g][[1]])){
        net$interactions[g][[1]][[i]]$probability <-
          net$interactions[g][[1]][[i]]$probability / c.sum
      }
    }
  }
  return(net)
}

##' Find the best regulatory input given the current regulatory function.
##'
##'
##' Try to add regulators to the "most recent" regulatory set. Stop only when
##' there is no improvement any longer.
##' The choice of the gene to add is done greedily.
##'
##' @title Greedily add new gene to regulatory function.
##' @param ts.multi A list of timeseries.
##' @param net A network.
##' @param g A gene.
##' @param partial Whether only the most recently added gene should have its
##'   parameter optimized. Defaults to FALSE.
##' @param verbal Whether messages to indicate progress should be shown.
##'   Defaults to TRUE.
##' @return A network.
innovateGene <- function(ts.multi, net, g,
                         partial=FALSE, verbal=TRUE){
  n <- nrow(ts.multi[[1]])
  iG <- net$interactions[[g]]
  L <- length(iG)
  iG <- iG[[L]]
  iG.orig <- iG
  iG.best <- iG

  ## This is what we have before doing anything
  losses <- rep(Inf, length(net$genes))
  loss.best <- LossGeneTotal(ts.multi, net, g)
  if (iG$input == 0){
    poss.add <- 1:n
  } else {
    poss.add <- (1:n)[-iG$input]
  }
  ## Try every possible addition of a regulator, and pick greedily
  for (reg in poss.add){
    if(iG.orig$input == 0){
      iG$input <- reg
    } else {
      iG$input <- c(iG.orig$input, reg)
    }
    net$interactions[[g]][[L]]$input <- iG$input
    k <- length(iG$input)
    if(partial){
      iG$a[k] <- bestPartialA(ts.multi, net, g, L, k)
    } else {
      iG$a <- bestA(ts.multi, net, g, L)
    }
    ## No point evaluating the loss function if we already know that the
    ## new gene doesn't add anything
    if(tail(iG$a, 1) != 0){
      net$interactions[[g]][[L]]$a <- iG$a
      func <- aToFunctionGene(net, g, L)
      net$interactions[[g]][[L]]$func <- func
      loss <- LossGeneTotal(ts.multi, net, g)
      ## Update if new best regulator set so far
      if (loss < loss.best){
        loss.best <- loss
        iG.best <- iG
        iG.best$func <- func
      }
    }
  }
  net$interactions[[g]][[L]] <- iG.best
  return(net)
}

##' Greedily add regulatory functions to explain a gene's time progression.
##'
##'
##' Genes are added until no more improvement is made or the upper limit is
##' reached. In biological networks, a common assumption for the an upper bound
##' on the number of inputs is in the range 3-5.
##'
##' @title Add new genes to regulatory set until no more improvement is made.
##' @param ts.multi A list of timeseries.
##' @param net A network.
##' @param g A gene.
##' @param partial Whether optimization of threshold parameter a is done only
##'   for latest gene. Defaults to FALSE.
##' @param verbal Whether progress messages should be displayed.
##' @param thresh By how much the improvement has to be to count. This has to be
##'   greater than zero as otherwise new genes will always be added to the brim.
##' @param k.max How many regulatory genes to allow maximally. To avoid
##'   excessive computation, this should be set to a small constant. Defaults to
##'   3.
##' @return A network.
innovateGeneUntilSaturated <- function(ts.multi, net, g,
                                       partial=FALSE, verbal=TRUE,
                                       thresh=0.1, k.max=3){
  ## Here we simply try adding new genes to the regulatory set until there is no
  ## benefit any longer.
  net.old <- net
  net.new <- net

  loss.new <- 0
  loss.old <- 1

  if(verbal){
    cat("Currently working on Gene", g, "\n")
  }
  k <- 1
  while ((k <= k.max) && (loss.new < loss.old - thresh)){
    net.old <- net.new
    loss.old <- LossGeneTotal(ts.multi, net.old, g)

    net.new <- innovateGene(ts.multi, net.new, g, partial)
    net.new <- optimizeC(ts.multi, net.new, g,
                         prune=FALSE)
    a.current <- tail(net.new$interactions[[g]], 1)[[1]]$a
    loss.new <- LossGeneTotal(ts.multi, net.new, g)
    k <- k + 1
  }
  return(net.new)
}

##' Add a new regulatory function complete with inputs.
##'
##'
##' This starts by taking some of the probability off the previous inferred
##' function and tries to find optimal genes and probabilities for the new
##' function.
##'
##' @title Add regulatory function to a gene.
##' @param ts.multi A list of timeseries.
##' @param net A network.
##' @param g A gene.
##' @param partial Whether optimization of threshold parameter a is done only
##'   for latest gene. Defaults to FALSE.
##' @param verbal Whether progress messages should be displayed.
##' @param thresh By how much the improvement has to be to count. This has to be
##'   greater than zero as otherwise new genes will always be added to the brim.
##' @param k.max How many regulatory genes to allow maximally. To avoid
##'   excessive computation, this should be set to a small constant. Defaults to
##'   3.
##' @param new.max An integer indicating how many new regulatory sets can be
##'   generated at once.
##' @return A network.
innovateGeneFunction <- function(ts.multi, net, g,
                                 partial=FALSE, verbal=TRUE,
                                 thresh=0.1, k.max=3,
                                 new.max=1){
  ## Tries to introduce a new regulatory set for a regulated gene.
  ## Not sure at all if this works yet.
  iter <- 1
  while(iter <= new.max){
    iG <- net$interactions[[g]]
    L <- length(iG)
    iG <- iG[[L]]
    iG$input <- 0
    iG$a <- 0
    iG$func <- c(-1, -1)
    net$interactions[[g]][[L]]$probability <- iG$probability/2
    iG$probability <- iG$probability/2
    net$interactions[[g]][[L+1]] <- iG
    net <- innovateGeneUntilSaturated(ts.multi, net, g,
                                      verbal=verbal, partial=partial,
                                      thresh=thresh, k.max=k.max)
    net <- optimizeC(ts.multi, net, g, prune=TRUE)
    iter <- iter + 1
  }
  return(net)
}
