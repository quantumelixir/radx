# find number of multi-indices preceding the given one
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
# around a point upto the dth order
# along -> matrix of multi-indices (one per row)
# S     -> restrict to a subspace spanned by rows of S
radxeval <- function(f, point, d=FALSE, along=FALSE, S=FALSE) {

    if (identical(d, FALSE) && identical(along, FALSE)) {
        stop('At least one of "d" or "along" must be specified.')
    }

    if (identical(d, FALSE)) {
        d = sum(along[1,])
    }

	# function from f: R^N -> R^M
	N <- length(point)
	M <- length(do.call(f, as.list(point)))

	# S is the NxN identity matrix by default
	if (identical(S, FALSE)) {
		S <- matrix(0, nrow=N, ncol=N)
		for (i in 1:N) 
			S[i, i] <- 1
	}
    Nr <- dim(S)[1]

    # If S is NrxN then J is PxNr
    J <- genMultiIndices(Nr, d)
    P <- dim(J)[1] # choose(Nr + d - 1, d)

	if (!identical(along, FALSE)) {
		# construct a partial interpolation matrix
		gammamatrix <- matrix(0, nrow=dim(along)[1], ncol=P)
		for (i in seq.int(dim(along)[1])) {
			gammamatrix[i,] <- gammaRow(along[i,])
		}
	}
	else {
        # construct the full interpolation matrix
		gammamatrix <- gammaMatrix(Nr, d)
	}

    # J:PxNr S:NrxN => rays:PxN
	rays <- J %*% S

	# Initialize the N P-directional d-order UTP variables
	x <- matrix(0, nrow=N, ncol=1 + P*d)
	x[,1] <- point
    indices <- 1 + d * seq(P) - (d - 1)
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

    derivatives <- matrix(nrow=P, ncol=M)

    # Propagate the UTPs!
    fp <- do.call(f, l)

    if (M > 1) {
        # for each scalar in the vector output
        for (m in seq.int(M)) {
            # pick the dth order derivatives
            dorder  <- fp[[m]]$coeff[1 + d * seq(P)]

            # interpolate to find the derivatives
            derivatives[, m] <- gammamatrix %*% dorder
        }
    }
    else {
            dorder  <- fp$coeff[1 + d * seq(P)]
            derivatives[,1] <- gammamatrix %*% dorder
    }

	return(derivatives)
}
