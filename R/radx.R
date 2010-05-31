# 1. Correctness
#   a. Compare with Finite Difference Methods
#   b. Richardson's Extrapolation
# 2. Optimization
#   a. Computer answer to each case explicity instead of delegating
#   b. Write code in C

radx <- function (val, type="variable", ord=1) {
    if (ord < 0 )
        stop("Error: radx objects must have ord >= 0")

    ad <- list()
    if (type == "variable") {
        if (ord == 0)
            ad$coeff <- c(val)
        else
            ad$coeff <- c(val, 1, rep(0, ord - 1))
    }
    else #constant
        ad$coeff <- c(val, rep(0, ord))

    ad$ord <- ord
    class(ad) <- "radx"

    return(ad)
}

radx_from <- function (coeffs) {
    ad <- list()
    ad$coeff <- c(coeffs)
    ad$ord <- length(ad$coeff) - 1
    class(ad) <- "radx"

    return(ad)
}

radx_empty <- function (ord) {
    if (ord < 0 )
        stop("Error: radx objects must have ord >= 0")

    return(radx_from(rep(0, ord + 1)))
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
    invisible(x)
}

