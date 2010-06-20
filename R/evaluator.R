# find number of multi-indices precedeing the given one
Nbefore <- function(J, n, v, cur=0) {
	if (n <= 1 || J[1] == v) # none precedes J in this case
		return(cur)

	for (i in v:(J[1]+1)) {
		cur <- cur + choose((n - 1) + (v - i) - 1, (v - i))
	}
	return(Nbefore(J[-1], n - 1, v - J[1], cur))
}

# each row of m (one multi-index) is converted to a position index
# returning a vector of positional indices (as many as rows in m)
convert2pos <- function(mat) {
	M <- dim(mat)[1]
	N <- dim(mat)[2]
    V <- sum(mat[1,])

    pos = rep(0, M) 
    for (m in seq(M))
		pos[m] <- 1 + Nbefore(mat[m,], N, V)

    return(pos)
}

# Compute derivatives of the function f
# around the point p upto the dth order
# only -> matrix of multi-indices (one per row)
# S    -> restrict to a subspace spanned by columns of S
radxeval <- function(f, p, d, only=FALSE, S=FALSE) {

	# function from f: R^N -> R^M
	N <- length(p)
	M <- length(do.call(f, as.list(p)))
	P <- choose(N + d - 1, d)

	# S is the NxN identity matrix by default
	if (identical(S, FALSE)) {
		S <- matrix(0, nrow=N, ncol=N)
		for (i in 1:N) 
			S[i, i] <- 1
	}

	# FIXME: BUGGY!!
	if (!identical(only, FALSE)) {
		# compute only specified derivatives
		J <- only

		# construct the interpolation matrix
		gammamatrix <- matrix(0, nrow=dim(only)[1], ncol=P)
		for (i in seq.int(dim(only)[1])) {
			gammamatrix[i,] <- gammaRow(only[i,])
		}

		indices <- 1 + d * convert2pos(only) - (d - 1)
	}
	else {
		# compute full derivative tensor
		J <- genMultiIndices(N, d)

		# construct the interpolation matrix
		gammamatrix <- gammaMatrix(N, d)

		indices <- 1 + d * seq(P) - (d - 1)
	}

	rays <- J %*% S

	# Initialize the N P-directional d-order UTP variables
	x <- matrix(0, nrow=N, ncol=1 + P*d)
	x[,1] <- p

	count <- 1
	for (i in indices) {
		x[,i] <- rays[count,]
		count <- count + 1
	}

	# pack the UTP variables into a list
	l <- list()
	for (i in seq.int(N)) {
		l[[i]] <- radx_from(x[i,], ndirs=P)
	}

	# Propagate the UTPs!
	fp <- do.call(f, l)
	print(fp)

	# pick the highest order derivatives
	dorder  <- fp$coeff[1 + d * seq(P)]

	# interpolate to find the derivatives
	derivatives <- gammamatrix %*% dorder

	return(derivatives)
}
