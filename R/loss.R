##' Compute loss of the network on the timeseries.
##'
##' This computes the complete loss of the network on a timeseries.
##' If desired, the prior penalty can be included.
##'
##' @param ts.multi List of timeseries.
##' @param net A PBN.
##' @param prior If TRUE, loss for the prior is included. Defaults to TRUE.
##' @return A float.
##' @seealso \code{\link{LossGeneTotal}} \code{\link{LossGeneTimeseries}} \code{\link{LossPriorGene}}
##' @examples
##' net <- createNetwork(inputProbabilities=1, 5, 2, "homogeneous")
##' ts.multi <- simulateNetwork(net, 10, 50) # Creates 50 timeseries à 10 points
##' loss <- Loss(ts.multi, net)
Loss <- function(ts.multi, net, prior=TRUE){
  ## The desired loss function we have here.
  loss <- 0
  loss.data <- 0
  loss.prior <- 0

  LossDF <- function(timeseries, net){
    loss.all <- t(sapply(1:nrow(timeseries),
                         function(g) LossGene(timeseries, net, g)))
    return(loss.all)
  }

  loss.measurement <- lapply(ts.multi, LossDF, net)
  loss.data <- Reduce("+", loss.measurement)

  loss <- sum(loss.data)
  if (prior){
    loss.prior <- LossPrior(net)
    loss <- loss + loss.prior
  }
  return(loss)
}


##' Compute the prior penalty for the network.
##'
##'
##' Can take either a BoolNet network or a TPBN and computes the penalty on the
##' network for its complexity.
##' The penalty is given by having to encode the number of inputs, the set of
##' actual inputs, and the Boolean function on these inputs.
##'
##' @param net A network, either threshold PBN or BoolNet net.
##' @param g The gene to be penalized.
##' @return A float.
##' @examples
##' net <- createNetwork(inputProbabilities=1, 5, 2, "homogeneous")
##' loss <- LossPriorGene(net, 1)
LossPriorGene <- function(net, g){
  n <- length(net$genes)
  interactions <- net$interactions[[g]]
  c.genes <- sapply(interactions, function(x) x$probability)

  ## If it's a threshold PBN, the a parameter is used. Otherwise the complexity
  ## of the parent functions are used.
  if(!is.null(interactions[[1]]$a)){
    A <- sapply(interactions, function(x) sum(x$a != 0))
    size.cost <- sum(log(1+A))
    choice.cost <- sum(log(choose(n, A)))
    function.cost <- sum(log(2^A))
  } else {
    parent.size <- sapply(interactions, function(x) length(x$input))
    size.cost <- sum(log(1+parent.size))
    choice.cost <- sum(log(choose(n, parent.size)))
    function.cost <- sum(log(2 ^ (2 ^ parent.size)))
  }
  return(size.cost + choice.cost + function.cost)
}

##' Computes the loss for a gene over a list of timeseries.
##'
##'
##' Makes use of \code{\link{LossGene}}, which computes the loss for only a
##' single timeseries.
##'
##' @param ts.multi A list of timeseries.
##' @param net A network.
##' @param g A gene.
##' @return A float.
##' @examples
##' net <- createNetwork(inputProbabilities=1, 5, 2, "homogeneous")
##' ts.multi <- simulateNetwork(net, 10, 50) # Creates 50 timeseries à 10 points
##' loss <- LossGeneTimeseries(ts.multi, net, 1)
LossGeneTimeseries <- function(ts.multi, net, g){
  losses <- sapply(ts.multi, function(ts) LossGene(ts, net, g))
  return(sum(losses))
}


##' Compute the total loss for a gene.
##'
##'
##' This includes both the predictions on the timeseries as well as the penalty
##' for its complexity.
##'
##' @param ts.multi A list of timeseries.
##' @param net A network.
##' @param g A gene.
##' @return A float.
##' @examples
##' net <- createNetwork(inputProbabilities=1, 5, 2, "homogeneous")
##' ts.multi <- simulateNetwork(net, 10, 50) # Creates 50 timeseries à 10 points
##' loss <- LossGeneTotal(ts.multi, net, 1)
LossGeneTotal <- function(ts.multi, net, g){
  return(LossGeneTimeseries(ts.multi, net, g) + LossPriorGene(net, g))
}


## These are all mostly just helper functions.

LossPrior <- function(net){
  n <- length(net$genes)
  loss.all <- t(sapply(1:n,
                       function(g) LossPriorGene(net, g)))
  return(sum(loss.all))
}

LossGene <- function(timeseries, net, g){
  gene.parents <- lapply(net$interactions[[g]], function(x) x$input)
  gene.functions <- lapply(net$interactions[[g]], function(x) x$func)
  gene.interactionProbabilities <- lapply(net$interactions[[g]],
                                          function(x) x$probability)
  p.flip <- net$flipProbability
  p.perturbation <- 1 - (1 - p.flip)^length(net$genes)
  p.correct <- 0
  loss <- 0
  for (t in 2:ncol(timeseries)){
    ## The new value of the gene might be due to a perturbation.
    if (timeseries[g, t] != timeseries[g, t-1]){
      p.correct <- p.correct + p.flip * p.perturbation
    } else {
      p.correct <- p.correct + (1 - p.flip) * p.perturbation
    }

    ## For each possible regulating function, determine whether it predicts
    ## the correct value for the gene
    for (i in 1:length(gene.parents)){
      parents <- gene.parents[[i]]
      funct <- gene.functions[[i]]
      c.genes <- gene.interactionProbabilities[[i]]
      p.correct <- p.correct + LossGene1(timeseries, t, g, i, net,
                                         parents, funct, c.genes)
    }
    loss <- loss - log(p.correct)

    ## Don't forget to reset the probability after each time step.
    p.correct <- 0
  }
  return(loss)
}

LossGene1 <- function(timeseries, t, g, i, net, parents, funct, c.genes){
  p.flip <- net$flipProbability
  p.perturbation <- 1 - (1 - p.flip)^length(net$genes)
  p.correct <- 0

  ## Compute value predicted by the i-th regulatory function for gene g.
  ## This may be of the form numeric(0), but that's just fine.
  parents <- na.omit(parents)
  parents.value <- timeseries[parents, t-1]
  gene.function.index <- unbinary(paste(parents.value, collapse="")) + 1
  gene.prediction <- funct[as.integer(gene.function.index)]
  if(gene.prediction == -1){
    gene.prediction <- timeseries[g, t-1]
  }

  ## If the prediction does match, its probability counts toward being correct.
  if (gene.prediction == timeseries[g, t]){
    p.correct <- ((1 - p.perturbation) * c.genes)
  }
  return(p.correct)
}
