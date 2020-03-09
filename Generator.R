require(mvtnorm)

SAVE_NODES = TRUE
SAVE_EDGES = TRUE
DATA_PREFIX = "A"
PATH = "/home/dgribel/src/ConstrainedSBM/data/"

save_samples <- function(data, omega, m, label) {
	# w_in  <- omega[1,1]
	# w_out <- omega[1,2]
	n <- nrow(data)
	d <- ncol(data)
	fileDesc <- paste(DATA_PREFIX, toString(m), toString(w_in), toString(w_out), sep = "-")
	fileDesc <- gsub('\\.', '', fileDesc)
	outputFile <- paste(PATH, fileDesc, sep = "")
    outputFile <- paste(outputFile, seed, sep = "-")
	outputFile <- paste(outputFile, "data", sep = ".")
	labelsFile <- paste(PATH, fileDesc, sep = "")
	labelsFile <- paste(labelsFile, seed, sep = "-")
    labelsFile <- paste(labelsFile, "label", sep = ".")
	# Write data to file
	# write.table(paste(toString(n), toString(d)), file = outputFile, sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
	# write.table(data, file = outputFile, sep = " ", row.names = FALSE, col.names = FALSE, append = FALSE)
	# Write labels to file
	write.table(label, file = labelsFile, sep = " ", row.names = FALSE, col.names = FALSE)
}

save_edges <- function(df, omega) {
	# w_in  <- omega[1,1]
	# w_out <- omega[1,2]
	fileDesc = paste(DATA_PREFIX, toString(m), toString(w_in), toString(w_out), sep = "-")
	fileDesc <- gsub('\\.', '', fileDesc)
	linksFile = paste(PATH, fileDesc, sep = "")
    linksFile = paste(linksFile, seed, sep = "-")
	linksFile = paste(linksFile, "link", sep = ".")
	cat(linksFile, "\n")
	# Write must-links to file
	write.table(df, file = linksFile, sep = " ", row.names = FALSE, col.names = FALSE)
}

generate_data <- function() {
	
	data <- matrix(nrow = n, ncol = d)

	means <- matrix(nrow = m, ncol = d)

	Sigma <- diag(1.25, d)

	# prb <- c(0.3, 0.2, 0.1, 0.4)

	# # Randomly choose the component that generates a data point 
	# label <- sample(seq(from = 1, to = m, by = 1), prob = prb, size = n, replace = TRUE)

	for(c in 1:m) {
		means[c,] <- rnorm(d, 0, 1)
	}

	################################# GENERATE SAMPLES

	for(i in 1:n) {
		data[i,] <- rmvnorm(1, means[label[i], ], Sigma)
	}

	################################# GENERATE EDGES


	# Omega = matrix( c(.040, .008,
	# 				  .008, .040),
	# 				nrow = m, ncol = m, byrow = TRUE )

	# Omega = matrix( c(.110, .010, .020,
	# 				  .010, .120, .015,
	# 				  .020, .015, .100),
	# 				nrow = m, ncol = m, byrow = TRUE )

	# Omega = matrix( c(.130, .010, .020, .015,
	# 				  .010, .140, .015, .013,
	# 				  .020, .015, .135, .011,
	# 				  .015, .013, .011, .120),
	# 				nrow = m, ncol = m, byrow = TRUE )

	edges <- data.frame()
	degree <- integer(n)

	for(i in 1:n) {
		y = 2*rpois(1, 0.5*Omega[label[i], label[i]])
		if(y > 0) {
			edges <- rbind(edges, c(i, i, y))
			degree[i] <- degree[i] + y
		}
	}

	for(i in 1:(n-1)) {
		for(j in (i+1):n) {
			y <- rpois(1, Omega[label[i], label[j]])
			if(y > 0) {
				edges <- rbind(edges, c(i, j, y))
				degree[i] <- degree[i] + y
				degree[j] <- degree[j] + y
			}
		}
	}

	# cat(Omega, "\n")
	save_samples(data, Omega, m, label)
	save_edges(edges, Omega)
}

args = commandArgs(trailingOnly = TRUE)

seed <- strtoi(args[1])

set.seed(seed)

n <- 100

m <- 4

# Number of features (dimensions)
d <- 2

diagonal = c(0.65, 0.7, 0.75, 0.8)
offDiagonal = rep(0.1, length(diagonal))

prb <- c(0.25, 0.25, 0.25, 0.25)

# Randomly choose the component that generates a data point 
label <- sample(seq(from = 1, to = m, by = 1), prob = prb, size = n, replace = TRUE)

# for(r in 1:(length(diagonal))) {
# 	w_in <- diagonal[r]
# 	w_out <- offDiagonal[r]
# 	Omega <- matrix(NA, ncol = m, nrow = m)
# 	Omega[lower.tri(Omega)] <- w_out
# 	Omega[upper.tri(Omega)] <- t(Omega)[upper.tri(t(Omega))]
# 	diag(Omega) <- w_in
# 	generate_data()
# }

w_in    <- 0.2
w_out   <- 0.06
eps_in  <- 0.1
eps_out <- 1.0

Omega <- matrix(NA, ncol = m, nrow = m)

for(r in 1:m) {
	Omega[r, r] <- w_in + runif(1, -eps_in*w_in, eps_in*w_in)
}

for(r in 1:(m-1)) {
	for(s in (r+1):m) {
		Omega[r, s] <- w_out + runif(1, -eps_out*w_out, eps_out*w_out)
		Omega[s, r] <- Omega[r, s]
	}
}

# Omega = matrix(   c(0.33848, 0.10403, 0.0845942, 0.0001,
# 					0.10403, 0.20196, 0.112381, 0.171528,
# 					0.0845942, 0.112381, 0.406302, 0.0001,
# 					0.0001, 0.171528, 0.0001, 0.36109),
# 					nrow = m, ncol = m, byrow = TRUE )

generate_data()