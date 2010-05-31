## Create a new algebra for "radx" objects
## Redefine basic arithmetic: +, -, *, /

"+.radx" <- function (x, y=NULL) {
    #unary plus
    if (is.null(y)) {
        return(x)
    }

    if (class(x) != "radx")
        return(y + x)

    coeff <- x$coeff
    if (class(y) == "radx") {
        if (x$ord != y$ord)
            stop("Error: '+' applied to radx objects of differing orders")
        coeff <- coeff + y$coeff
    }
    else { #numeric
        coeff[1] <- coeff[1] + y
    }

    return(radx_from(coeff))
}

"-.radx" <- function (x, y=NULL) {
    # unary negation
    if (is.null(y)) {
        return(radx_from(-x$coeff))
    }

    if (class(x) != "radx")
        return(-(y-x))

    coeff <- x$coeff
    if (class(y) == "radx") {
        if (x$ord != y$ord)
            stop("Error: '-' applied to radx objects of differing orders")

        coeff <- coeff - y$coeff
    }
    else { #numeric
        coeff[1] <- coeff[1] - y
    }

    return(radx_from(coeff))
}

"*.radx" <- function (x, y) {
    if (class(x) != "radx") {
        coeff <- x * y$coeff
    }
    else if (class(y) == "radx") {
        if (x$ord != y$ord)
            stop("Error: '*' applied to radx objects of differing orders")

        coeff <- rep(0, x$ord + 1)
        coeff[1] <- x$coeff[1] * y$coeff[1]
        if(x$ord)
            for (k in seq(x$ord)) {
                i <- 0
                while (i <= k) {
                    coeff[k + 1] <- coeff[k + 1] + x$coeff[i + 1] * y$coeff[k - i + 1] / (base::factorial(i) * base::factorial(k - i))
                    i <- i + 1
                }
                coeff[k + 1] <- coeff[k + 1] * base::factorial(k)
            }
    }
    else { #numeric
        coeff <- y * x$coeff
    }

    return(radx_from(coeff))
}

"/.radx" <- function (x, y) {
    if (class(x) != "radx") {
        # replace with specific code for constant / y
        return(radx_from(c(x, rep(0, y$ord))) / y)
    }
    else if (class(y) == "radx") {
        if (x$ord != y$ord)
            stop("Error: '/' applied to radx objects of differing orders")

        coeff <- rep(0, x$ord + 1)
        coeff[1] <- x$coeff[1] / y$coeff[1]
        for (k in seq(x$ord)) {
            i <- 0
            while (i < k) {
                coeff[k + 1] <- coeff[k + 1] + coeff[i + 1] * y$coeff[k - i + 1] /(base::factorial(i) * base::factorial(k - i))
                i <- i + 1
            }
            coeff[k + 1] <- base::factorial(k) * (x$coeff[k + 1] - coeff[k + 1]) / y$coeff[1]
        }
    }
    else { #numeric
        coeff <- x$coeff / y
    }

    return(radx_from(coeff))
}

