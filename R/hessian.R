hess <- function (val, type="variable", vars=1, wdir=1) {
    if (vars < 1)
        stop("Error: hess objects must have vars >= 1")

    if (wdir < 1 || wdir > vars)
        stop("Error: wdir must be >= 1 and <= vars")

    hess <- list()
    hess$grad <- c(val, rep(0, vars))
    hess$hess <- matrix(0, nrow=vars, ncol=vars)

    if (type == "variable")
        hess$grad[wdir + 1] <- 1

    class(hess) <- "hess"

    return(hess)
}

hess_from <- function (g, h) {
    hess <- list()
    hess$grad<- c(g)
    hess$hess<- h
    class(hess) <- "hess"
    return(hess)
}

hess_empty <- function (vars=1) {
    if (vars < 1)
        stop("Error: hess objects must have vars >= 1")

    return(hess_from(rep(0, vars + 1), matrix(0, vars, vars)))
}

print.hess <- function (x) {
    cat("Value      :", x$grad[1], "\n")

    if (length(x$grad) > 1) {
        cat("Gradient:")
        print.default(format(x$grad[-1]), print.gap = 2, quote = FALSE)
        cat("Hessian:\n")
        print(x$hess)
    }
    else
        cat("Gradient: NA\nHessian: NA\n")

    invisible(x)
}

# f   = u + v
# f'  = u' + v'
# f'` = u'` + v'`
"+.hess" <- function (x, y=NULL) {
    #unary plus
    if (is.null(y)) {
        return(x)
    }

    if (class(x) != "hess")
        return(y + x)
    else if (class(y) == "hess") {
        if (length(x$grad) != length(y$grad)) 
            stop("Error: '+' applied to hess objects of differing number of variables")
        grad <- x$grad + y$grad
        hess <- x$hess + y$hess
    }
    else {
        grad <- x$grad
        grad[1] <- grad[1] + y
        hess <- x$hess
    }

    return(hess_from(grad, hess))
}

# f   = u - v
# f'  = u' - v'
# f'` = u'` - v'`
"-.hess" <- function (x, y=NULL) {
    #unary minus
    if (is.null(y)) {
        return(x)
    }

    if (class(x) != "hess")
        return(-(y - x))
    else if (class(y) == "hess") {
        if (length(x$grad) != length(y$grad))
            stop("Error: '-' applied to hess objects of differing number of variables")
        grad <- x$grad - y$grad
        hess <- x$hess - y$hess
    }
    else {
        grad <- x$grad
        grad[1] <- grad[1] - y
        hess <- x$hess
    }

    return(hess_from(grad, hess))
}

# f   = u * v
# f'  = u'v + uv'
# f'` = u'`v + u`v' + u'v` + uv'`
"*.hess" <- function (x, y) {
    if (class(x) != "hess") {
        return(y * x)
    }
    else if (class(y) == "hess") {
        if (length(x$grad) != length(y$grad))
            stop("Error: '*' applied to hess objects of differing number of variables")

        l <- length(x$grad) - 1
        grad <- rep(0, length(x$grad))
        grad[1] <- x$grad[1] * y$grad[1]
        for (i in seq(l)) {
            grad[i + 1] <- x$grad[i + 1] * y$grad[1] + x$grad[1] * y$grad[i + 1]
        }

        hess <- matrix(nrow=l, ncol=l)
        for (i in seq(l))
            for (j in seq(l))
                hess[i, j] <- x$hess[i, j] * y$grad[1] + x$grad[i + 1] * y$grad[j + 1] + x$grad[j + 1] * y$grad[i + 1] + x$grad[1] * y$hess[i, j]
    }
    else {
        grad <- y * x$grad
        hess <- y * x$hess
    }

    return(hess_from(grad, hess))
}

# f   = u / v
# f'  = (u'v - uv')/v^2
# f'` = ((u'`v + u'v` - u`v' - uv'`)(v^2) - (2vv`)(u'v - uv')) / (v^4)
"/.hess" <- function (x, y) {
    if (class(x) != "hess") {
        return(hess_from(c(x, rep(0, length(y$grad) - 1)), matrix(0, nrow=length(y$grad)-1, ncol=length(y$grad)-1)) / y)
    }
    else if (class(y) == "hess") {
        if (length(x$grad) != length(y$grad))
            stop("Error: '/' applied to hess objects of differing number of variables")

        l <- length(x$grad) - 1
        grad <- rep(0, length(x$grad))
        grad[1] <- x$grad[1] / y$grad[1]
        for (i in seq(l)) {
            grad[i + 1] <- (x$grad[i + 1] * y$grad[1] - x$grad[1] * y$grad[i + 1]) / (y$grad[1] * y$grad[1])
        }

        hess <- matrix(nrow=l, ncol=l)
        for (i in seq(l))
            for (j in seq(l))
                hess[i, j] <- ((x$hess[i,j]*y$grad[1] + x$grad[i+1]*y$grad[j+1] - x$grad[j+1]*y$grad[i+1] - x$grad[1]*y$hess[i,j])*(y$grad[1]^2) - (2*y$grad[1]*y$grad[j+1])*(x$grad[i+1]*y$grad[1] - x$grad[1]*y$grad[i+1])) / (y$grad[1] ^ 4)
    }
    else {
        grad <- x$grad / y
        hess <- x$hess / y
    }

    return(hess_from(grad, hess))
}
