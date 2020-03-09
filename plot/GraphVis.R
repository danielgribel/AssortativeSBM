library(GGally)
library(network)
library(sna)
library(ggplot2)

seed = 1001
set.seed(seed)

INSTANCE <- "Network561"

FOLDER <- "/home/dgribel/src/ConstrainedSBM/data/"

PATH_EDGES <- paste(FOLDER, INSTANCE, ".link", sep = "")
PATH_LABEL <- paste(FOLDER, INSTANCE, ".label", sep = "")
# PATH_TRUTH <- paste(FOLDER, INSTANCE, ".truth", sep = "")

edges <- read.table(PATH_EDGES, header = FALSE)
label <- read.table(PATH_LABEL, header = FALSE)
# truth <- read.table(PATH_TRUTH, header = FALSE)

n <- nrow(label)

A <- matrix(0, ncol = n, nrow = n)

for (i in 1:nrow(edges)) {
    a <- edges[i, 1]
    b <- edges[i, 2]
    w <- edges[i, 3]
    A[a, b] <- w
    A[b, a] <- w
}

net = network(A, directed = FALSE)

# vertex names
network.vertex.names(net) = letters[1:n]

x <- paste(as.vector(t(label)))

x <- replace(x, x == "1", "blue")
x <- replace(x, x == "2", "red")
x <- replace(x, x == "3", "green")
x <- replace(x, x == "4", "orange")
x <- replace(x, x == "5", "yellow")
x <- replace(x, x == "6", "gray")
x <- replace(x, x == "7", "black")

# pdf("Cortex-DC-CSBM2.pdf")

g1 <- ggnet2(net, label = 1:n, node.size = 6, node.color = x, edge.size = 1, edge.color = "grey")

plot(g1)

# dev.off()