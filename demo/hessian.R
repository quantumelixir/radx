library(radx)

f <- function(x,y,z) {
    x*(x*(x+y) / z)
}

N <- 3
d <- 2

J <- genMultiIndices(N, d)

S <- matrix(0, nrow=N, ncol=N)
for (i in 1:N) 
    S[i, i] <- 1

rays <- J %*% S
P <- dim(rays)[1]

point <- c(3,3,3)

x <- radx(point[1], ord=d, ndirs=P, ray=rays[,1])
y <- radx(point[2], ord=d, ndirs=P, ray=rays[,2])
z <- radx(point[3], ord=d, ndirs=P, ray=rays[,3])

a <- f(x, y, z)

dth_order  <- a$coeff[seq.int(1 + a$ord, 1 + a$ord*a$ndirs, a$ord)]

gammamatrix <-gammaMatrix(N, d)

partials <- gammamatrix %*% dth_order

H <- matrix(0, nrow=N, ncol=N)
k <- 1
for (i in 1:N) {
    for (j in i:N) {
        H[i, j] <- H[j, i] <- partials[k]
        k <- k + 1
    }
}

cat('Hessian of the function: f(x, y, z) = x*(x*(x+y)/z) at (3,3,3) is:\n')
H
