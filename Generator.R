require(mvtnorm)

DATA_PREFIX = "C"
PATH = "/home/dgribel/src/ConstrainedSBM/data/"

saveLabels <- function(label, omega) {
	fileDesc <- paste(DATA_PREFIX, toString(m), toString(wIn), toString(wOut), sep = "-")
	fileDesc <- gsub('\\.', '', fileDesc)
	labelsFile <- paste(PATH, fileDesc, sep = "")
	labelsFile <- paste(labelsFile, seed, sep = "-")
	labelsFile <- paste(labelsFile, "label", sep = ".")
	write.table(label, file = labelsFile, sep = " ", row.names = FALSE, col.names = FALSE)
	cat('File created at:', labelsFile, '\n')
}

saveEdges <- function(network, omega) {
	fileDesc = paste(DATA_PREFIX, toString(m), toString(wIn), toString(wOut), sep = "-")
	fileDesc <- gsub('\\.', '', fileDesc)
	linksFile = paste(PATH, fileDesc, sep = "")
	linksFile = paste(linksFile, seed, sep = "-")
	linksFile = paste(linksFile, "link", sep = ".")
	write.table(network, file = linksFile, sep = " ", row.names = FALSE, col.names = FALSE)
	cat('File created at:', linksFile, '\n')
}

createNetwork <- function() {
	edges <- data.frame()

	# Sample A_ii from a Poisson distribution
	for(i in 1:n) {
		a_ii = 2*rpois(1, 0.5*Omega[label[i], label[i]])
		if(a_ii > 0) {
			edges <- rbind(edges, c(i, i, a_ii))
		}
	}

	# Sample A_ij, for i != j, from a Poisson distribution
	for(i in 1:(n-1)) {
		for(j in (i+1):n) {
			a_ij <- rpois(1, Omega[label[i], label[j]])
			if(a_ij > 0) {
				edges <- rbind(edges, c(i, j, a_ij))
			}
		}
	}
	saveLabels(label, Omega)
	saveEdges(edges, Omega)
}

seed <- 1
set.seed(seed)

# Number of samples
n <- 100

# Number of communities
m <- 4

# Randomly choose the community that generates a sample
label <- sample(seq(from = 1, to = m, by = 1), prob = rep(1/m, m), size = n, replace = TRUE)

wIn    <- 0.5
wOut   <- 0.2
epsIn  <- 0.1
epsOut <- 1.0

Omega <- matrix(NA, ncol = m, nrow = m)

# Define the values of \omega in the diagonal
for(r in 1:m) {
	Omega[r, r] <- wIn + runif(1, -epsIn*wIn, epsIn*wIn)
}

# Define the values of \omega in the off-diagonal
for(r in 1:(m-1)) {
	for(s in (r+1):m) {
		Omega[r, s] <- wOut + runif(1, -epsOut*wOut, epsOut*wOut)
		Omega[s, r] <- Omega[r, s]
	}
}

createNetwork()