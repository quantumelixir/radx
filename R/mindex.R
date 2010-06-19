# the global multi-index matrix with one mult-index per row
mat <- matrix()

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

# Generate multi-indices of given length (N) and absolute value (V)
genMultiIndices <- function(N, V) {
    mat <<- matrix(nrow=choose(N + V - 1, V), ncol=N)
    fillMultiIndices(1, 1, N, V)
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
gammaMultiIndex <- function(i, j, d) {
    if (length(i) != length(j))
        stop('Error: Lengths of the multi-index arguments are unequal')
    if (absolute(i) != absolute(j))
        stop('Error: Absolute values of the multi-index arguments are unequal')

    n <- length(i)
    d <- absolute(i)

    gij <- 0
    for (a in 1:d) {
        genMultiIndices(n, a)
        total <- dim(mat)[1]
        for (b in 1:total) {
            k <- mat[b,]
            gij <- gij + singleTerm(i, j, k, d)
        }
    }
    return(gij)
}

# perform testing
gtest <- function(N, V) {
    genMultiIndices(N, V)
    mindices <- dim(mat)[1]
    gammamatrix <- matrix(0, nrow=mindices, ncol=mindices)

    nonzero <- 0
    for (i in 1:mindices) {
    for (j in 1:mindices) {
        gammamatrix[i,j] <- gammaMultiIndex(mat[i,], mat[j,])
        if (abs(gammamatrix[i, j]) > 1e-6)
            nonzero <- nonzero + 1
    }
    }

    mval <- 0
    for (i in 1:mindices) {
        cur <- sum(abs(gammamatrix[i,]))
        if (mval < cur)
            mval = cur
    }

    cat('N =                ', N, '\n')
    cat('D =                ', V, '\n')
    cat('Nonzero Entries   :', nonzero, '\n')
    cat('Max_i (Sum_j g_ij):', mval, '\n')
    print(gammamatrix)
}
