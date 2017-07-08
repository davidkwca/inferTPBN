require(BoolNet)
require(entropy)
require(foreach)
require(doParallel)
require(igraph)
require(compositions)
require(Rsolnp)
##' Helper function for testing network inference.
##'
##' Automatically generate a network, generate timeseries data from it, and do inference.
##' @title
##' @param n
##' @param k
##' @param p
##' @param num.timepoints
##' @param num.experiments
##' @param topology
##' @param gamma
##' @param n.cores
##' @param seed
##' @param partial
##' @param verbal
##' @return
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

setupPBN <- function(n=20, k=5, p=0.01,
                     inputProbabilities=c(0.5, 0.4, 0.1),
                     topology="homogeneous",
                     num.timepoints=10, num.experiments=20,
                     seed=111,
                     verbal=FALSE){
  ## Create a network with the desired number n of genes, having k inputs on average.
  ## The topology determines whether the network is homogeneous, scale-free, or
  set.seed(seed)
  net.true <- createNetwork(inputProbabilities=inputProbabilities,
                            n=n, k=k,
                            topology=topology,
                            noIrrelevantGenes=TRUE,
                            readableFunctions=TRUE)
  net.true <- transformNetwork(net.true)

  ts.multi <- simulateNetwork(net.true, num.timepoints, num.experiments)
  return(list(net.true=net.true, ts.multi=ts.multi))
}

## for (n in c(5, 10, 25, 50)){
##   for (k in c(3, 5, 10)){
##     for (e in c(20, 50, 100)){
##       tryCatch({
##         FullRun(n=n, k=k, num.experiments=e, n.cores=detectCores())
##       }, error=function(cond){sprintf("Error for n=%d, k=%d, e=%d", n, k, e)})
##     }
##   }
## }

## for (n in c(5, 10, 25, 50, 100)){
##   FullRun(n=n, k=5, num.experiments=100, n.cores=detectCores()-1)
## }

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

ev = FALSE
if(ev){
  for (n in c(10, 25, 50)){
    l <- load.nets(n=n, num.experiments=50, k=5)
    net.true <- l$net.true
    net.inferred <- l$net.inferred
    s <- EvaluateStatistics(net.true, net.inferred)
    cat(sprintf("Size %d. Recall: %.2f,  Precision: %.2f, F-score: %.2f, Structural Accuracy: %.2f", n, s$recall, s$precision, s$F, s$str.acc), "\n")
  }
}
## for (g in 1:n){
##   s <- EvaluateStatisticsGene(net.true, net.inferred, g)
##   cat(sprintf("Gene %d. Recall: %.2f,  Precision: %.2f, F-score: %.2f",
##               g, s$recall, s$precision, s$F), "\n")
## }
