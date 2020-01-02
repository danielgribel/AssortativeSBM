require(mvtnorm)

SAVE_NODES = TRUE
SAVE_EDGES = TRUE
DATA_PREFIX = "Data"
PATH = "/home/dgribel/src/ConstrainedSBM/data/"

seed = 1000
set.seed(seed)

save_samples <- function(m, data, label) {
	n <- nrow(data)
	d <- ncol(data)
	fileDesc = paste(DATA_PREFIX, toString(m), toString(d), toString(n), sep = "-")
	outputFile = paste(PATH, fileDesc, sep = "")
    outputFile = paste(outputFile, seed, sep = "-")
	outputFile = paste(outputFile, "data", sep = ".")
	labelsFile = paste(PATH, fileDesc, sep = "")
	labelsFile = paste(labelsFile, seed, sep = "-")
    labelsFile = paste(labelsFile, "label", sep = ".")
	# Write data to file
	write.table(paste(toString(n), toString(d)), file = outputFile, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(data, file = outputFile, sep = " ", row.names = FALSE, col.names = FALSE, append=TRUE)
	# Write labels to file
	write.table(label, file = labelsFile, sep = " ", row.names = FALSE, col.names = FALSE)
}

save_edges <- function(n, m, df) {
	fileDesc = paste(DATA_PREFIX, toString(m), toString(d), toString(n), sep = "-")
	linksfile = paste(PATH, fileDesc, sep = "")
    linksfile = paste(linksfile, seed, sep = "-")
	linksfile = paste(linksfile, "link", sep = ".")
	# Write must-links to file
	# write.table(paste(toString(n), toString(nrow(df))), file = linksfile, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
	write.table(df, file = linksfile, sep = " ", row.names = FALSE, col.names = FALSE, append=TRUE)
}

# Number of data points
n <- 100

# Number of Gaussian components (clusters)
m <- 2

# Number of features (dimensions)
d <- 8

data <- matrix(nrow = n, ncol = d)

means <- matrix(nrow = m, ncol = d)

Sigma <- diag(1.25, d)

# Randomly choose the component that generates a data point 
label <- sample(seq(from = 1, to = m, by = 1), size = n, replace = TRUE)

for(c in 1:m) {
    means[c,] <- rnorm(d, 0, 1)
}

################################# GENERATE SAMPLES

for(i in 1:n) {
    data[i,] <- rmvnorm(1, means[label[i], ], Sigma)
}

################################# GENERATE EDGES

# Assortative networks, n = 200
# omega_in = 0.015
# omega_out = 0.001

# Dis-Assortative networks, n = 200
# omega_in = 0.001
# omega_out = 0.025

# Assortative networks, n = 400
# omega_in = 0.0075
# omega_out = 5e-04

# Dis-Assortative networks, n = 400
# omega_in = 5e-04
# omega_out = 0.0125

Omega = matrix( c(.075, .002, .002, .075), nrow = m, ncol = m, byrow = TRUE )

edges <- data.frame()

for(i in 1:n) {
	y = 2*rpois(1, 0.5*Omega[label[i], label[i]])
	if(y > 0) {
		edges <- rbind(edges, c(i, i, y))
	}
}

for(i in 1:(n-1)) {
	for(j in (i+1):n) {
        y <- rpois(1, Omega[label[i], label[j]])
		if(y > 0) {
			edges <- rbind(edges, c(i, j, y))
		}
	}
}

save_samples(m, data, label)
save_edges(n, m, edges)