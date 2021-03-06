# 1. Correctness
#   a. Compare with Finite Difference Methods
#   b. Richardson's Extrapolation
# 2. Optimization
#   a. Computer answer to each case explicity instead of delegating
#   b. Write code in C

radx <- function (val, ord=1, ndirs=1, ray=TRUE) {
    if (ord < 0 || ndirs < 1)
        stop("Error: radx objects must have ord >= 0, dirs >= 1")

    ad <- list()
    ad$ord <- ord
    ad$ndirs <- ndirs

    if (ord == 0)
        ad$coeff <- c(val)
    else
    {
        ad$coeff <- c(val, rep(rep(0, ord), ndirs))
        if (!isTRUE(ray)) {
            ad$coeff[seq.int(2, length(ad$coeff) - ad$ord + 1, ad$ord)] <- ray
        }
        else if (ad$ord > 0) {
            ad$coeff[2] <- 1
        }
    }

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

# generate list of factorials starting from 1:
# 1!, ..., n! (inclusive of n)
genfact <- function(n) {
    ret <- rep(1, n)
    if (n > 1)
        for (i in 2:n) {
            ret[i] <- ret[i - 1] * i
        }
    return(ret)
}

print.radx <- function (x) {
    cat("Value      :", x$coeff[1], "\n")

    tcoeff <- x$coeff[-1]
    ftlist <- genfact(x$ord)

    if (x$ord > 0) {
        cat("Coefficients:")
        print.default(format(tcoeff), print.gap = 2, quote = FALSE)
    }
    else
        cat("Coefficients: NA\n")

    if (x$ord > 0) {
        cat("Derivatives :")
        print.default(format(tcoeff * ftlist), print.gap = 2, quote = FALSE)
    }
    else
        cat("Derivatives : NA\n")

    cat("Order       : ", x$ord, "\n", sep="")
    cat("Directions  : ", x$ndirs, "\n", sep="")
    invisible(x)
}

# overload concatenation for radx objects
c.radx <- function(...) {
    return(list(...))
}

# get the n-th order derivative of the m-th output function
getDeriv <- function(x, n, m=1) {
    if (class(x) != "radx")
        stop("Error: getDeriv accepts only radx objects")

    if (n > x$ord)
        stop("Error: Requesting an order higher than computed")
    if (n < 0)
        stop("Error: Bad value for derivative order")

    if (m > 1) {
        coeff <- x$coeff[n + 1,m]
    }
    else {
        coeff <- x$coeff[n + 1]
    }

    return(coeff * factorial(n))
}
