# 1. Correctness
#   a. Compare with Finite Difference Methods
#   b. Richardson's Extrapolation
# 2. Optimization
#   a. Computer answer to each case explicity instead of delegating
#   b. Write code in C

radx <- function (val, type="variable", ord=1, ndirs=1, wdir=1) {
    if (ord < 0 || ndirs < 1)
        stop("Error: radx objects must have ord >= 0, dirs >= 1")

    if (wdir > ndirs)
        stop("Error: wdir greater than ndirs")

    ad <- list()
    if (type == "variable") {
        if (ord == 0)
            ad$coeff <- c(val)
        else
            ad$coeff <- c(val, rep(rep(0, ord), ndirs))
            ad$coeff[2 + ord * (wdir - 1)] <- 1
    }
    else if (type == "constant")
        ad$coeff <- c(val, rep(rep(0, ord), ndirs))
    else
        stop('Bad type: Choose between "variable" and "constant"')

    ad$ord <- ord
    ad$ndirs <- ndirs

    class(ad) <- "radx"

    return(ad)
}

radx_from <- function (coeffs, ndirs=1) {
    ad <- list()
    ad$ndirs <- ndirs
    ad$coeff <- c(coeffs)
    ad$ord <- (length(ad$coeff) - 1)/ndirs
    class(ad) <- "radx"

    return(ad)
}

radx_empty <- function (ord, ndirs=1) {
    if (ord < 0 || ndirs < 1)
        stop("Error: radx objects must have ord >= 0, dirs >= 1")

    return(radx_from(c(0, rep(rep(0, ord), ndirs)), ndirs))
}

print.radx <- function (x) {
    cat("Value      :", x$coeff[1], "\n")

    if (x$ord > 0) {
        cat("Derivatives:")
        print.default(format(x$coeff[-1]), print.gap = 2, quote = FALSE)
    }
    else
        cat("Derivatives: NA\n")

    cat("Order      : ", x$ord, "\n", sep="")
    cat("Directions : ", x$ndirs, "\n", sep="")
    invisible(x)
}

