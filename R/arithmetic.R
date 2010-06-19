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
        if (x$ord != y$ord || x$ndirs != y$ndirs)
            stop("Error: '+' applied to radx objects of differing orders and/or directions")
        coeff <- coeff + y$coeff
    }
    else { #numeric
        coeff[1] <- coeff[1] + y
    }

    return(radx_from(coeff, ndirs=x$ndirs))
}

"-.radx" <- function (x, y=NULL) {
    # unary negation
    if (is.null(y)) {
        return(radx_from(-x$coeff, ndirs=x$ndirs))
    }

    if (class(x) != "radx")
        return(-(y-x))

    coeff <- x$coeff
    if (class(y) == "radx") {
        if (x$ord != y$ord || x$ndirs != y$ndirs)
            stop("Error: '-' applied to radx objects of differing orders and/or directions")

        coeff <- coeff - y$coeff
    }
    else { #numeric
        coeff[1] <- coeff[1] - y
    }

    return(radx_from(coeff, ndirs=x$ndirs))
}

"*.radx" <- function (x, y) {
    if (class(x) != "radx") {
        coeff <- x * y$coeff
        return(radx_from(coeff, ndirs=y$ndirs))
    }
    else if (class(y) == "radx") {
        if (x$ord != y$ord || x$ndirs != y$ndirs)
            stop("Error: '*' applied to radx objects of differing orders and/or directions")

        coeff <- rep(0, length(x$coeff))
        coeff[1] <- x$coeff[1] * y$coeff[1]
        if(x$ord)
            for (p in seq(x$ndirs))
                for (d in seq(x$ord)) {
                    s <- (p - 1) * x$ord

                    coeff[s + d + 1] <- x$coeff[1] * y$coeff[s + d + 1]
                    i <- 1
                    while (i < d) {
                        coeff[s + d + 1] <- coeff[s + d + 1] + x$coeff[s + i + 1] * y$coeff[s + d - i + 1]
                        i <- i + 1
                    }
                    coeff[s + d + 1] <- coeff[s + d + 1] + x$coeff[s + d + 1] * y$coeff[1]
                }
    }
    else { #numeric
        coeff <- y * x$coeff
    }

    return(radx_from(coeff, ndirs=x$ndirs))
}

"/.radx" <- function (x, y) {
    if (class(x) != "radx") {
        # replace with specific code for constant / y
        return(radx_from(c(x, rep(0, length(y$coeff) - 1)), ndirs=y$ndirs) / y)
    }
    else if (class(y) == "radx") {
        if (x$ord != y$ord || x$ndirs != y$ndirs)
            stop("Error: '/' applied to radx objects of differing orders and/or directions")

        coeff <- rep(0, length(x$coeff))
        coeff[1] <- x$coeff[1] / y$coeff[1]
        if(x$ord)
            for (p in seq(x$ndirs))
                for (d in seq(x$ord)) {
                    s <- (p - 1) * x$ord

                    coeff[s + d + 1] <- coeff[1] * y$coeff[s + d + 1]
                    i <- 1
                    while (i < d) {
                        coeff[s + d + 1] <- coeff[s + d + 1] + coeff[s + i + 1] * y$coeff[s + d - i + 1]
                        i <- i + 1
                    }

                    coeff[s + d + 1] <- (x$coeff[s + d + 1] - coeff[s + d + 1]) / y$coeff[1]
                }
    }
    else { #numeric
        coeff <- x$coeff / y
    }

    return(radx_from(coeff, ndirs=x$ndirs))
}

