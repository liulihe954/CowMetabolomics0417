library(qgraph)
data(big5)
data(big5groups)
str(big5)
str(big5groups)

set.seed(1) 
adj <- matrix(sample(0:1,10^2,TRUE,prob=c(0.8,0.2)),nrow=10,ncol=10) 
Q <- qgraph(adj) 
centrality(Q)


set.seed(1) 
x <- cor(matrix(rnorm(25), nrow = 5)) 
colors <- c("red", "red", "blue", "blue", "white") 

# colored qgraph plot 
qg <- qgraph(x, colors = colors) 
# randomly assing motifs to colors (notice that white nodes stay white) 
makeBW(qg) 
# associate a motif only to one of the colors 
makeBW(qg, colorlist = c("blue")) 
# define an order, which allows to choose motifs 
makeBW(qg, colorlist = c("blue", "red"),plot = TRUE)
makeBW(qg, colorlist = c("red", "blue"),plot = TRUE)
