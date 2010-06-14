grad <- function (val, type="variable", vars=1, wdir=1) {
    if (vars < 1)
        stop("Error: grad objects must have vars >= 1")

    if (wdir < 1 || wdir > vars)
        stop("Error: wdir must be >= 1 and <= vars")

    grad <- list()
    grad$coeff <- c(val, rep(0, vars))

    if (type == "variable")
        grad$coeff[wdir + 1] <- 1

    class(grad) <- "grad"

    return(grad)
}

grad_from <- function (coeffs) {
    grad <- list()
    grad$coeff <- c(coeffs)
    class(grad) <- "grad"
    return(grad)
}

grad_empty <- function (vars=1) {
    if (vars < 1)
        stop("Error: grad objects must have vars >= 1")

    return(grad_from(rep(0, vars + 1)))
}

print.grad <- function (x) {
    cat("Value      :", x$coeff[1], "\n")

    if (length(x$coeff) > 1) {
        cat("Gradient:")
        print.default(format(x$coeff[-1]), print.gap = 2, quote = FALSE)
    }
    else
        cat("Gradient: NA\n")

    invisible(x)
}

# f  = u + v
# f' = u' + v'
"+.grad" <- function (x, y=NULL) {
    #unary plus
    if (is.null(y)) {
        return(x)
    }

    if (class(x) != "grad")
        return(y + x)
    else if (class(y) == "grad") {
        if (length(x$coeff) != length(y$coeff)) 
            stop("Error: '+' applied to grad objects of differing number of variables")
            coeff <- x$coeff + y$coeff
    }
    else {
        coeff <- x$coeff
        coeff[1] <- coeff[1] + y
    }

    return(grad_from(coeff))
}

# f  = u - v
# f' = u' - v'
"-.grad" <- function (x, y=NULL) {
    #unary minus
    if (is.null(y)) {
        return(x)
    }

    if (class(x) != "grad")
        return(y + x)
    else if (class(y) == "grad") {
        if (length(x$coeff) != length(y$coeff)) 
            stop("Error: '-' applied to grad objects of differing number of variables")
            coeff <- x$coeff - y$coeff
    }
    else {
        coeff <- x$coeff
        coeff[1] <- coeff[1] - y
    }

    return(grad_from(coeff))
}

# f  = u * v
# f' = u'v + uv'
"*.grad" <- function (x, y) {
    if (class(x) != "grad") {
        return(y * x)
    }
    else if (class(y) == "grad") {
        if (length(x$coeff) != length(y$coeff))
            stop("Error: '*' applied to grad objects of differing number of variables")

        l <- length(x$coeff) - 1
        coeff <- rep(0, length(x$coeff))
        coeff[1] <- x$coeff[1] * y$coeff[1]
        for (i in seq(l)) {
            coeff[i + 1] <- x$coeff[i + 1] * y$coeff[1] + x$coeff[1] * y$coeff[i + 1]
        }
    }
    else {
        coeff <- y * x$coeff
    }

    return(grad_from(coeff))
}

# f  = u / v
# f' = (u'v - uv')/v^2
"/.grad" <- function (x, y) {
    if (class(x) != "grad") {
        return(grad_from(c(x, rep(0, length(y$coeff) - 1))) / y)
    }
    else if (class(y) == "grad") {
        if (length(x$coeff) != length(y$coeff))
            stop("Error: '/' applied to grad objects of differing number of variables")

        l <- length(x$coeff) - 1
        coeff <- rep(0, length(x$coeff))
        coeff[1] <- x$coeff[1] / y$coeff[1]
        for (i in seq(l)) {
            coeff[i + 1] <- (x$coeff[i + 1] * y$coeff[1] - x$coeff[1] * y$coeff[i + 1]) / (y$coeff[1] * y$coeff[1])
        }
    }
    else {
        coeff <- x$coeff / y
    }

    return(grad_from(coeff))
}
