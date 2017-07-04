plotStats <- function(n=c(5, 10, 25, 50, 100), k=c(5), p=c(0.01),
                     num.timepoints=c(10), num.experiments=c(1),
                     topology="homogeneous"){
  N <- length(n)
  K <- length(k)
  P <- length(p)
  T <- length(num.timepoints)
  E <- length(num.experiments)
  R <- length(topology)


  time.complete <- array(dim=c(N, K, P, T, E, R))
  time.bestfit <- array(dim=c(N, K, P, T, E, R))
  loss.complete <- array(dim=c(N, K, P, T, E, R))
  loss.bestfit <- array(dim=c(N, K, P, T, E, R))
  loss.true <- array(dim=c(N, K, P, T, E, R))
  for (i in 1:N){
    for (j in 1:K){
      for (l in 1:P){
        for (t in 1:T){
          for (e in 1:E){
            for (r in 1:R){
              dir <- sprintf("points-%d_experiments-%d_topology-%s/",
                             num.timepoints[t], num.experiments[e], topology[r])
              file.name <- sprintf("%sn-%d_k-%d_p-%.2f", dir, n[i], k[j], p[l])
              load(file=file.name)
              time.complete[i, j, l, t, e, r] <- inferred.list$time.complete['elapsed']
              loss.complete[i, j, l, t, e, r] <-
                Loss(setup.list$ts.multi,
                     inferred.list$net.inferred,
                     prior=FALSE)
              loss.true[i, j, l, t, e, r] <-
                Loss(setup.list$ts.multi,
                     setup.list$net.true,
                     prior=FALSE)
            }
          }
        }
      }
    }
  }
  ## By SIZE.
  if (N > 1){
    times <- time.complete[, 1, 1, 1, 1, 1]
    png(filename=sprintf("times_by_size.png"))
    plot(x=n, y=times,
         xlab="network size",
         ylab="inference time")
    title("Time")
    dev.off()

    loss.inferred <- loss.complete[, 1, 1, 1, 1, 1]
    loss.true <- loss.true[, 1, 1, 1, 1, 1]
    ## loss.bestfit <- loss.bestfit[, 1, 1, 1, 1, 1]

    yLim <- c(0, max(c(loss.true, loss.inferred, loss.bestfit)))
    png(filename=sprintf("losses_by_size.png"))
    plot(x=n, y=loss.true, col="red", ylim=yLim, type="l",
         xlab="network size",
         ylab="loss")
    lines(x=n, y=loss.inferred, col="blue", type="l")
    legend("topleft", c("true", "inferred", "bestfit"),
           lty=1, col=c("red", "blue", "green"))
    title("Losses")
    dev.off()
    png(filename=sprintf("relative_loss_by_size.png"))
    plot(x=n, y=(loss.inferred-loss.true)/loss.true)
    dev.off()
  }
  ## By EXPERIMENTS.
  if (E > 1){
    times <- time.complete[1, 1, 1, 1, , 1]
    png(filename=sprintf("times_by_experiments.png"))
    plot(x=num.experiments, y=times,
         xlab="num experiments",
         ylab="inference time")
    title("Time")
    dev.off()

    loss.true <- loss.true[1, 1, 1, 1, , 1]
    loss.inferred <- loss.complete[1, 1, 1, 1, , 1]
    loss.bestfit <- loss.bestfit[1, 1, 1, 1, , 1]
    yLim <- c(0, max(c(loss.true, loss.inferred, loss.bestfit)))

    png(filename=sprintf("losses_by_experiments.png"))
    plot(x=num.experiments, y=loss.true, col="red", ylim=yLim, type="l",
         xlab="num experiments",
         ylab="loss")
    lines(x=num.experiments, y=loss.inferred, col="blue", type="l")
    lines(x=num.experiments, y=loss.bestfit, col="green", type="l")
    legend("topleft", c("true", "inferred", "bestfit"),
           lty=1, col=c("red", "blue", "green"))
    title("Losses")
    dev.off()

    png(filename=sprintf("relative_loss_by_experiments.png"))
    plot(x=num.experiments, y=(loss.inferred-loss.true)/loss.true)
    dev.off()
  }
}

EvaluateStatisticsGene <- function(net.true, net, g){
  returns <- list()
  n <- length(net.true$genes)
  interactions.true <- net.true$interactions[[g]]
  interactions <- net$interactions[[g]]
  inputs.true <- lapply(interactions.true, function(x) x$input)
  inputs <- lapply(interactions, function(x) x$input)
  inputs.true.flat <- unlist(inputs.true)
  inputs.flat <- unlist(inputs)

  noninputs.true.flat <- setdiff(1:10, inputs.true.flat)
  noninputs.flat <- setdiff(1:10, inputs.flat)

  num.interaction.true <- length(unique(inputs.true.flat)) # TP + FN
  num.interaction.found <- length(unique(inputs.flat)) # TP + FP
  num.interaction.intersection <- length(unique(intersect(inputs.true.flat, inputs.flat)))

  num.noninteraction.true <- length(unique(noninputs.true.flat))
  num.noninteraction.found <- length(unique(noninputs.flat))
  num.noninteraction.intersection <-
    length(unique(intersect(noninputs.true.flat, noninputs.flat)))

  recall <- num.interaction.intersection / num.interaction.true
  precision <- num.interaction.intersection / num.interaction.found

  returns$num.interaction.true <- num.interaction.true
  returns$num.interaction.found <- num.interaction.found
  returns$num.interaction.intersection <- num.interaction.intersection
  returns$num.noninteraction.true <- num.noninteraction.true
  returns$num.noninteraction.found <- num.noninteraction.found
  returns$num.noninteraction.intersection <- num.noninteraction.intersection

  returns$recall <- recall
  returns$precision <- precision
  if(all(precision == 0, recall == 0)){
    returns$F <- 0
  } else {
    returns$F <- 2 * (precision * recall)/(precision + recall)
  }
  return(returns)
}

EvaluateStatistics <- function(net.true, net){
  n <- length(net$genes)

  interaction.true.total <- 0
  interaction.found.total <- 0
  interaction.intersection.total <- 0

  noninteraction.true.total <- 0
  noninteraction.found.total <- 0
  noninteraction.intersection.total <- 0


  F.avg <- 0
  for (g in 1:n){
    stats <- EvaluateStatisticsGene(net.true, net, g)
    interaction.true.total <- interaction.true.total + stats$num.interaction.true
    interaction.found.total <- interaction.found.total + stats$num.interaction.found
    interaction.intersection.total <- interaction.intersection.total + stats$num.interaction.intersection

    noninteraction.true.total <- noninteraction.true.total + stats$num.noninteraction.true
    noninteraction.found.total <- noninteraction.found.total + stats$num.noninteraction.found
    noninteraction.intersection.total <- noninteraction.intersection.total + stats$num.noninteraction.intersection
  }

  recall.avg <- interaction.intersection.total / interaction.true.total
  precision.avg <- interaction.intersection.total / interaction.found.total
  F.avg <- 2 * (precision.avg * recall.avg)/(precision.avg + recall.avg)
  str.acc <- (interaction.intersection.total + noninteraction.intersection.total) / n^2

  return(list(recall=recall.avg, precision=precision.avg, F=F.avg, str.acc=str.acc))
}

PrintLosses <- function(ts.multi, nets){
  cat("Loss by genes:", "\n")
  for (g in 1:n){
    cat("Gene:", g, "\n")
    i <- 0
    for (net in nets){
      i <- i+1
      cat(sprintf("Likelihood of net %d:", i),
          LossGeneTimeseries(ts.multi, net, g), "| ")
      cat(sprintf("Prior of net %d:", i),
          LossPriorGene(net, g),
          "\n")
    }
  }
}

CompareInteractions <- function(nets){
  ## Useless so far.
  for (g in 1:n){
    input.list[[g]] <- list()
    input.list[[g]][[1]] <- net.true$interactions[[g]][[1]]$input
    input.list[[g]][[2]] <- net.inferred$interactions[[g]][[1]]$input
    input.list[[g]][[3]] <- net.inferred.partial$interactions[[g]][[1]]$input
    input.list[[g]][[4]] <- net.bestfit$interactions[[g]][[1]]$input
  }
}

PBNtoGraph <- function(net, fname=""){
  n <- length(net$genes)
  edges <- c()
  weight <- c()
  for (g in 1:n){
    parents.g <- lapply(net$interactions[[g]], function(x) x$input)

    prob.g <- rank(sapply(net$interactions[[g]], function(x) x$probability))
    for (i in 1:length(prob.g)){
      parents <- parents.g[[i]]
      edges <- c(edges, Reduce(c, sapply(parents, function(x) c(x, g))))
      weight <- c(weight, rep(prob.g[i], length(parents)))
    }
  }
  G <- graph(edges)
  E(G)$weight <- weight

  q <- sort(unique(floor(rank(weight))))
  for (i in 1:length(q))
    E(G)$color[floor(rank(E(G)$weight)) == q[i]] <- i
  if(fname != ""){
    pdf(file=fname)
    plot(g.40, vertex.size=4, edge.arrow.size=0.5, edge.arrow.width=0.3, vertex.size2=1, vertex.label.cex=0.4)
    dev.off()
  }
  return(G)
}
