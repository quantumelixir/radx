# UTP Algebra for Elementary Functions:
# exp, log, a ^ x, x ^ a, x ^ x, sin, cos

"exp.radx" <- function(x) {
    if(class(x) == "radx") {
        coeff <- rep(0, length(x$coeff))
        coeff[1] <- exp(x$coeff[1])
        for (p in seq(x$ndirs))
            for (d in seq(x$ord)) {
                s <- (p - 1) * x$ord

                i <- 1
                while (i < d) {
                    coeff[s + d + 1] <- coeff[s + d + 1] + i * x$coeff[s + i + 1] * coeff[s + d - i + 1] / (base::factorial(i) * base::factorial(d - i))
                    i <- i + 1
                }
                coeff[s + d + 1] <- coeff[s + d + 1] + d * x$coeff[s + d + 1] * coeff[1]
                coeff[s + d + 1] <- coeff[s + d + 1] / d
            }
    }
    else { #numeric
        return(exp(x))
    }

    return(radx_from(coeff, ndirs=x$ndirs))
}

# natural logarithm (base e)
"log.radx" <- function(x) {
    if(class(x) == "radx") {
        coeff <- rep(0, length(x$coeff))
        coeff[1] <- log(x$coeff[1])
        for (p in seq(x$ndirs))
            for (d in seq(x$ord)) {
                s <- (p - 1) * x$ord

                i <- 1
                while (i < d) {
                    coeff[s + d + 1] <- coeff[s + d + 1] + i * coeff[s + i + 1] * x$coeff[s + d - i + 1]
                    i <- i + 1
                }

                coeff[s + d + 1] <- (x$coeff[s + d + 1] - coeff[s + d + 1] / d) / x$coeff[1]
            }
    }
    else { #numeric
        return(log(x))
    }

    return(radx_from(coeff, ndirs=x$ndirs))
}

# a ^ x and x ^ a
"^.radx" <- function(x, y) {
    if(class(x) == "radx" && class(y) == "radx") { # x ^ x
        return(exp(y * log(x)))
    }
    else if(class(x) == "radx") {
        coeff <- rep(0, length(x$coeff))
        coeff[1] <- x$coeff[1] ^ y
        for (p in seq(x$ndirs))
            for (d in seq(x$ord)) {
                s <- (p - 1) * x$ord

                i <- 1
                while (i < d) {
                    coeff[s + d + 1] <- coeff[s + d + 1] + ((y + 1) ^ i / d - 1) * x$coeff[s + i + 1] * coeff[s + d - i + 1]
                    i <- i + 1
                }
                coeff[s + d + 1] <- coeff[s + d + 1] + ((y + 1) ^ d / d - 1) * x$coeff[s + d + 1] * coeff[1]

                coeff[s + d + 1] <- coeff[s + d + 1] / x$coeff[1]
            }
    }
    else if(class(y) == "radx") {
        return(exp(y * log(x)))
    }
    else { #numeric
        return(x ^ y)
    }

    return(radx_from(coeff, ndirs=x$ndirs))
}

"sin.radx" <- function(x) {
    if(class(x) == "radx") {
        s_coeff <- rep(0, length(x$coeff))
        c_coeff <- rep(0, length(x$coeff))

        s_coeff[1] <- sin(x$coeff[1])
        c_coeff[1] <- cos(x$coeff[1])
        for (p in seq(x$ndirs))
            for (d in seq(x$ord)) {
                s <- (p - 1) * x$ord

                i <- 1
                while (i < d) {
                    s_coeff[s + d + 1] <- s_coeff[s + d + 1] + i * x$coeff[s + i + 1] * c_coeff[s + d - i + 1]
                    c_coeff[s + d + 1] <- c_coeff[s + d + 1] + i * x$coeff[s + i + 1] * s_coeff[s + d - i + 1]
                    i <- i + 1
                }
                s_coeff[s + d + 1] <- s_coeff[s + d + 1] + d * x$coeff[s + d + 1] * c_coeff[1]
                c_coeff[s + d + 1] <- c_coeff[s + d + 1] + d * x$coeff[s + d + 1] * s_coeff[1]

                s_coeff[s + d + 1] <- s_coeff[s + d + 1] / d
                c_coeff[s + d + 1] <- -1 * c_coeff[s + d + 1] / d
            }
    }
    else { #numeric
        return(sin(x))
    }

    return(radx_from(s_coeff, ndirs=x$ndirs))
}

"cos.radx" <- function(x) {
    if(class(x) == "radx") {
        s_coeff <- rep(0, length(x$coeff))
        c_coeff <- rep(0, length(x$coeff))

        s_coeff[1] <- sin(x$coeff[1])
        c_coeff[1] <- cos(x$coeff[1])
        for (p in seq(x$ndirs))
            for (d in seq(x$ord)) {
                s <- (p - 1) * x$ord

                i <- 1
                while (i < d) {
                    s_coeff[s + d + 1] <- s_coeff[s + d + 1] + i * x$coeff[s + i + 1] * c_coeff[s + d - i + 1]
                    c_coeff[s + d + 1] <- c_coeff[s + d + 1] + i * x$coeff[s + i + 1] * s_coeff[s + d - i + 1]
                    i <- i + 1
                }
                s_coeff[s + d + 1] <- s_coeff[s + d + 1] + d * x$coeff[s + i + 1] * c_coeff[1]
                c_coeff[s + d + 1] <- c_coeff[s + d + 1] + d * x$coeff[s + i + 1] * s_coeff[1]

                s_coeff[s + d + 1] <- s_coeff[s + d + 1] / d
                c_coeff[s + d + 1] <- -1 * c_coeff[s + d + 1] / d
            }
    }
    else { #numeric
        return(cos(x))
    }

    return(radx_from(c_coeff, ndirs=x$ndirs))
}

"tan.radx" <- function(x) {
    if(class(x) == "radx") {
        p_coeff <- rep(0, length(x$coeff))
        w_coeff <- rep(0, length(x$coeff))

        p_coeff[1] <- tan(x$coeff[1])
        w_coeff[1] <- 1 + tan(x$coeff[1])^2
        for (p in seq(x$ndirs))
            for (d in seq(x$ord)) {
                s <- (p - 1) * x$ord

                i <- 1
                while (i < d) {
                    p_coeff[s + d + 1] <- p_coeff[s + d + 1] + i * x$coeff[s + i + 1] * w_coeff[s + d - i + 1]
                    w_coeff[s + d + 1] <- w_coeff[s + d + 1] + i * p_coeff[s + i + 1] * p_coeff[s + d - i + 1]
                    i <- i + 1
                }
                p_coeff[s + d + 1] <- (p_coeff[s + d + 1] + d * x$coeff[s + d + 1] * w_coeff[1]) / d
                w_coeff[s + d + 1] <- 2 * (w_coeff[s + d + 1] + d * p_coeff[s + d + 1] * p_coeff[1]) / d
            }
    }
    else { #numeric
        return(tan(x))
    }

    return(radx_from(p_coeff, ndirs=x$ndirs))
}

# BUGGY! FIXME: asin(x) computes correctly upto the second order, fails afterwards. Why?!
"asin.radx" <- function(x) {
    if(class(x) == "radx") {
        p_coeff <- rep(0, length(x$coeff))
        w_coeff <- rep(0, length(x$coeff))

        p_coeff[1] <- asin(x$coeff[1])
        w_coeff[1] <- sqrt(1 - x$coeff[1]^2)
        for (p in seq(x$ndirs))
            for (d in seq(x$ord)) {
                s <- (p - 1) * x$ord

                i <- 1
                while (i < d) {
                    p_coeff[s + d + 1] <- p_coeff[s + d + 1] + i * p_coeff[s + i + 1] * w_coeff[s + d - i + 1]
                    w_coeff[s + d + 1] <- w_coeff[s + d + 1] + i * p_coeff[s + i + 1] * x$coeff[s + d - i + 1]
                    i <- i + 1
                }
                p_coeff[s + d + 1] <- (d * x$coeff[s + d + 1] - p_coeff[s + d + 1])/(d * w_coeff[1])
                w_coeff[s + d + 1] <- -1 * (w_coeff[s + d + 1] + d * p_coeff[s + d + 1] * x$coeff[1]) / d
            }
    }
    else { #numeric
        return(asin(x))
    }

    return(radx_from(p_coeff, ndirs=x$ndirs))
}

"atan.radx" <- function(x) {
    if(class(x) == "radx") {
        p_coeff <- rep(0, length(x$coeff))
        w_coeff <- rep(0, length(x$coeff))

        p_coeff[1] <- atan(x$coeff[1])
        w_coeff[1] <- (1 + x$coeff[1]^2)
        for (p in seq(x$ndirs))
            for (d in seq(x$ord)) {
                s <- (p - 1) * x$ord

                i <- 1
                while (i < d) {
                    p_coeff[s + d + 1] <- p_coeff[s + d + 1] + i * p_coeff[s + i + 1] * w_coeff[s + d - i + 1]
                    w_coeff[s + d + 1] <- w_coeff[s + d + 1] + i * x$coeff[s + i + 1] * x$coeff[s + d - i + 1]
                    i <- i + 1
                }
                p_coeff[s + d + 1] <- (d * x$coeff[s + d + 1] - p_coeff[s + d + 1])/(d * w_coeff[1])
                w_coeff[s + d + 1] <- 2 * (w_coeff[s + d + 1] + d * x$coeff[s + d + 1] * x$coeff[1]) / d
            }
    }
    else { #numeric
        return(atan(x))
    }

    return(radx_from(p_coeff, ndirs=x$ndirs))
}

"acot" <- function(x) {
    return(-atan(x))
}
