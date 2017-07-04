##' Returns the loss of the network on the timeseries.
##'
##' Stuff.
##' @title Loss function for time series.
##' @param ts.multi List of timeseries.
##' @param net A PBN.
##' @param prior whether or not the prior should be included. Defaults to TRUE.
##' @return A float.
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


##' @title Penalty for single gene's regulatory set's complexity.
##' @param net A network, either threshold PBN or BoolNet net.
##' @param g The gene to be penalized.
##' @return A float.
LossPriorGene <- function(net, g){
  n <- length(net$genes)
  interactions <- net$interactions[[g]]
  c.genes <- sapply(interactions, function(x) x$probability)
  entropy <- entropy.empirical(c.genes)

  ## If it's threshold PBN, the a parameter is used. Otherwise the complexity
  ## of the parent functions are used.
  if(!is.null(interactions[[1]]$a)){
    A <- sapply(interactions, function(x) sum(x$a != 0))
    choice.cost <- sum(log(choose(n, A)))
    function.cost <- sum(log(2^A))
  } else {
    parent.size <- sapply(interactions, function(x) length(x$input))
    choice.cost <- sum(log(choose(n, parent.size)))
    function.cost <- sum(log(2 ^ (2 ^ parent.size)))
  }
  return(choice.cost + function.cost + entropy)
}


##' @title The loss on timeseries data for a single gene.
##' @param ts.multi A list of timeseries.
##' @param net A network.
##' @param g A gene.
##' @return A float.
LossGeneTimeseries <- function(ts.multi, net, g){
  losses <- sapply(ts.multi, function(ts) LossGene(ts, net, g))
  return(sum(losses))
}


##' @title The total loss of the network on a single gene.
##' @param ts.multi A list of timeseries.
##' @param net A network.
##' @param g A gene.
##' @return A float.
LossGeneTotal <- function(ts.multi, net, g){
  return(LossGeneTimeseries(ts.multi, net, g) + LossPriorGene(net, g))
}


## These are all mostly just helper functions.

LossPrior <- function(net){
  n <- length(net$genes)
  loss.all <- t(sapply(1:n,
                       function(g) LossPriorGene(net, g)))
  return(loss.all)
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
