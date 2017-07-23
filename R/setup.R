require(BoolNet)
require(entropy)
require(foreach)
require(doParallel)
require(igraph)
require(compositions)
require(Rsolnp)

##' Complete run of network generation and inference.
##'
##'
##' Automatically generate a network, generate timeseries data from it, and do inference.
##' Saves the inferred and true networks to a subdirectory.
##'
##' @param n Size of the network.
##' @param k The number of inputs per regulatory function for each gene, if homogeneous
##'   topology is used
##' @param p The probability of a perturbation.
##' @param num.timepoints The number of time points per timeseries generated.
##' @param num.experiments The number of timeseries to generate.
##' @param topology The topology to be used. Can be "homogeneous" or "scale_free".
##' @param gamma The exponent for the power law if topology = "scale_free".
##' @param n.cores The number of cores to use in the inference.
##' @param seed The random seed to use.
##' @param partial If TRUE, a network using partial optimization should be inferred. Defaults to FALSE.
##' @param verbal If TRUE, show progress as to which genes are currently being worked on.
##' @seealso \code{\link{setupPBN}} \code{\link{inferPBN}}
FullRun <- function(n=20, k=5, p=0.01,
                    num.timepoints=10, num.experiments=50,
                    topology="homogeneous", gamma=2.5,
                    n.cores=detectCores()-1,
                    seed=111,
                    partial=FALSE, verbal=FALSE){

  dir <- sprintf("points-%d_experiments-%d_topology-%s/",
                 num.timepoints, num.experiments, topology)
  dir.create(dir, showWarnings=FALSE)
  file.name <-sprintf("%sn-%d_k-%d_p-%.2f", dir, n, k, p)
  if(file.exists(file.name)){
    cat(sprintf("File %s already exists!", file.name), "\n")
    return()
  } else {
    setup.list <- setupPBN(n=n, k=k, p=p,
                           num.timepoints=num.timepoints,
                           num.experiments=num.experiments,
                           topology=topology,
                           seed=seed)
    net.true <- setup.list$net.true
    ts.multi <- setup.list$ts.multi

    inferred.list <- inferPBN(ts.multi,
                              n.cores=n.cores,
                              seed=seed,
                              partial=partial,
                              verbal=verbal)

    save(setup.list, inferred.list, n.cores,
         file=file.name)
  }
}
##' Generate a random TPBN and a timeseries simulated from it.
##'
##'
##' Creates a network with the desired characteristics and simulates a list of timeseries
##' from it which can be used to test the network inference algorithm against.
##'
##' @param n Size of the network.
##' @param k The number of inputs per regulatory function for each gene, if homogeneous
##'   topology is used.
##' @param p The probability of a perturbation.
##' @param inputProbabilities Probabilities for creating one, two, three, ... regulatory
##'   functions for a gene.
##' @param topology The topology to be used. Can be "homogeneous" or "scale_free".
##' @param gamma The exponent for the power law if topology = "scale_free".
##' @param num.timepoints The number of time points per timeseries generated.
##' @param num.experiments The number of timeseries to generate.
##' @param seed The random seed to use.
##' @return A list containing the generated network and a list of timeseries simulated
##'   from it.
##' @seealso \code{\link{transformNetwork}} \code{\link{simulateNetwork}}
##' @examples
##' setupPBN(n=5, k=3, topology="scale_free", gamma=2.5)
setupPBN <- function(n=20, k=5, p=0.01,
                     inputProbabilities=c(0.5, 0.4, 0.1),
                     topology="homogeneous", gamma=2.5,
                     num.timepoints=10, num.experiments=20,
                     seed=111,
                     verbal=FALSE){
  ## Create a network with the desired number n of genes, having k inputs on average.
  ## The topology determines whether the network is homogeneous, scale-free, or
  set.seed(seed)
  net.true <- createNetwork(inputProbabilities=inputProbabilities,
                            n=n, k=k,
                            topology=topology, gamma=gamma,
                            noIrrelevantGenes=TRUE,
                            readableFunctions=TRUE)
  net.true <- transformNetwork(net.true)

  ts.multi <- simulateNetwork(net.true, num.timepoints, num.experiments)
  return(list(net.true=net.true, ts.multi=ts.multi))
}

load.nets <- function(n=25, k=5,
                      num.experiments=20, num.timepoints=10,
                      p=0.01, topology="homogeneous"){
  dir <- sprintf("points-%d_experiments-%d_topology-%s/",
                 num.timepoints, num.experiments, topology)
  file.name <- sprintf("%sn-%d_k-%d_p-%.2f", dir, n, k, p)
  load(file.name)
  return(list(net.true=setup.list$net.true,
              net.inferred=inferred.list$net.inferred))
}
