# Generate multi-indices of given length (N) and absolute value (V)
genMultiIndices <- function(N, V) {
    mat <- matrix(nrow=choose(N + V - 1, V), ncol=N)

    # fillMultiIndices: Helper function for genMultiIndices
    # Fill choose(n + v - 1, v) multi-indices of length n
    # and sum v in the matrix mat starting from (i,j)
    fillMultiIndices <- function(i, j, n, v) {
        if (n == 1) {
            mat[i, j] <<- v
        }
        else if (v == 0) {
            mat[i, j:(j + n - 1)] <<- 0
        }
        else {
            rowOffset <- 0
            # the first element of each multi-index can be any of 0, 1, ..., v
            for (k in v:0) {
                times <- choose((n - 1) + (v - k) - 1, (v - k))
                mat[(i + rowOffset):(i + rowOffset + times - 1), j] <<- k
                fillMultiIndices(i + rowOffset, j + 1, n - 1, v - k)
                rowOffset <- rowOffset + times
            }
        }
    }

    fillMultiIndices(1, 1, N, V)
    return(mat)
}

# compute choose(z, k) where z and k are multi-indices
binomialMultiIndex <- function(z, k) {
    prod(choose(z, k))
}

# absolute value of multi-index objects
absolute <- function(z) {
    sum(z)
}

# compute one term of the multi-index summation
singleTerm <- function(i, j, k, d) {
    t1 <- (-1) ^ absolute(i - k)
    t2 <- binomialMultiIndex(i, k)
    t3 <- binomialMultiIndex(d * k / absolute(k), j)
    t4 <- (absolute(k) / d) ^ absolute(i)
    return(t1 * t2 * t3 * t4)
}

# compute the exact interpolation coefficients
gammaMultiIndex <- function(i, j) {
    if (length(i) != length(j))
        stop('Error: Lengths of the multi-index arguments are unequal')
    if (absolute(i) != absolute(j))
        stop('Error: Absolute values of the multi-index arguments are unequal')

    n <- length(i)
    d <- absolute(i)

    gij <- 0
    for (a in 1:d) {
        J <- genMultiIndices(n, a)
        total <- dim(J)[1]
        for (b in 1:total) {
            k <- J[b,]
            gij <- gij + singleTerm(i, j, k, d)
        }
    }
    return(gij)
}

# computes the interpolation coefficients for calculating
# the partial derivative corresponding to a multi-index i
gammaRow <- function(i) {
    N <- length(i)
    V <- absolute(i)
    P <- choose(N + V - 1, V)

    J <- genMultiIndices(N, V)
    gammarow <- matrix(nrow=1, ncol=P)
    for (j in 1:P) {
        gammarow[j] <- gammaMultiIndex(i, J[j,])
    }
    return(gammarow)
}

# gamma matrix stores interpolation coefficients for all
# (pure and mixed) derivatives of order V in N variables
gammaMatrix <- function(N, V) {
    J <- genMultiIndices(N, V)
    P <- choose(N + V - 1, V)
    gammamatrix <- matrix(0, nrow=P, ncol=P)
    for (i in 1:P) {
        for (j in 1:P) {
            gammamatrix[i,j] <- gammaMultiIndex(J[i,], J[j,])
        }
    }
    return(gammamatrix)
}
