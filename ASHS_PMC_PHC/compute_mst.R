#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(igraph)

adjlist=read.table(args[1]);
N=max(adjlist[,c(1,2)]);

adj=matrix(data=NA,nrow=N,ncol=N);
for (i in 1:dim(adjlist)[1]) {
  adj[adjlist[i,1],adjlist[i,2]] = 1.0 - adjlist[i,3];
}

g=graph_from_adjacency_matrix((adj + t(adj)),weighted='weight',diag=F,mode='undirected')

T=minimum.spanning.tree(g, algorithm='prim')

D=distances(g, algorithm='dijkstra')

# Compute the distances from D
root=which.min(colSums(D));

# Print the shortest paths
path=shortest_paths(T,from=root)
for (i in 1:N) {
  if(i == root)
    cat(i,"\n")
  else
    cat(rev(path$vpath[[i]]),"\n")
}
