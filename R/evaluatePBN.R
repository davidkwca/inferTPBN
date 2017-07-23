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

FscoreNet <- function(net.true, net.inferred){
  n <- length(net.true$genes)

  tp <- 0
  fp <- 0
  predicted <- 0
  relevant <- 0

  for (g in 1:n){
    inter <- net.inferred$interactions[[g]]
    input <- c()
    for (j in 1:length(inter)){
      input <- c(input, inter[[j]]$input)
    }
    input.true <- c()
    for (j in 1:length(net.true$interactions[[g]])){
      input.true <- c(input.true, net.true$interactions[[g]][[j]]$input)
    }
    predicted <- predicted + length(unique(input))
    relevant <- relevant + length(unique(input.true))
    tp <- tp + length(unique(intersect(input, input.true)))
  }
  precision <- tp / predicted
  recall <- tp / relevant
  f <- 2 * (precision * recall) / (precision + recall)
  return(list(precision=precision, recall=recall, f=f))
}

loadNetStats <- function(n, k, e=20, t=10, topology='homogeneous'){
  dir <- sprintf("points-%d_experiments-%d_topology-%s/",
                 t, e, topology)
  file.name <- sprintf("%sn-%d_k-%d_p-%.2f", dir, n, k, 0.01)
  load(file=file.name)

  net.true <- setup.list$net.true
  ## ts.multi <- setup.list$ts.multi

  ts.multi <- simulateNetwork(net.true, 10, 50)

  net.inferred <- inferred.list$net.inferred
  time.inferred <- inferred.list$time.complete

  loss.true <- Loss(ts.multi, net.true, prior=FALSE)
  loss.inferred <- Loss(ts.multi, net.inferred, prior=FALSE)

  f=FscoreNet(net.true, net.inferred)

  return(list(loss.true=loss.true,
              loss.inferred=loss.inferred,
              time.inferred=time.inferred,
              f=f))
}

plotStats <- function(ns, k, e=20, t=10, topology='homogeneous'){
  losses.true <- rep(0, length(ns))
  losses.inferred <- rep(0, length(ns))
  fs <- rep(0, length(ns))
  times <- rep(0, length(ns))
  for (i in 1:length(ns)){
    n <- ns[i]
    l <- loadNetStats(n, k, e, t, topology)
    losses.true[i] <- l$loss.true
    losses.inferred[i] <- l$loss.inferred
    times[i] <- l$time.inferred[3]
    fs[i] <- l$f
  }
  yLim <- c(min(losses.true), max(losses.inferred))
  png(sprintf("plots/losses-k-%d_e-%d_topology-%s.png", k, e, topology))
  plot(x=ns, y=losses.true, col="red", type="l", ylim=yLim,
       xlab="network size",
       ylab="loss")
  lines(x=ns, y=losses.inferred, col="blue", type="l")
  if (e == 20){
    legend("topleft", legend=c("True", "Inferred"),
       col=c("red", "blue"), lty=c(1, 1), cex=1.7)
  }
  dev.off()

  png(sprintf("plots/relative-k-%d_e-%d_topology-%s.png", k, e, topology))
  plot(x=ns, y=(losses.inferred-losses.true)/losses.true)
  dev.off()

  png(sprintf("plots/fscore-k-%d_e-%d_topology-%s.png", k, e, topology))
  plot(x=ns, y=fs)
  dev.off()

  png(sprintf("plots/times-k-%d_e-%d_topology-%s.png", k, e, topology))
  plot(x=ns, y=times, log="xy", type="l", col="black")
  dev.off()
}
