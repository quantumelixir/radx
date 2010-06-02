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
                coeff[s + d + 1] <- coeff[s + d + 1] + d * x$coeff[s + d + 1] * coeff[1] / base::factorial(d)

                coeff[s + d + 1] <- coeff[s + d + 1] * base::factorial(d - 1)
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
                    coeff[s + d + 1] <- coeff[s + d + 1] + i * coeff[s + i + 1] * x$coeff[s + d - i + 1] / (base::factorial(i) * base::factorial(d - i))
                    i <- i + 1
                }

                coeff[s + d + 1] <- (x$coeff[s + d + 1] - base::factorial(d - 1) * coeff[s + d + 1]) / x$coeff[1]
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
                    coeff[s + d + 1] <- coeff[s + d + 1] + ((y + 1) ^ i / d - 1) * x$coeff[s + i + 1] * coeff[s + d - i + 1] / (base::factorial(i) * base::factorial(d - i))
                    i <- i + 1
                }
                coeff[s + d + 1] <- coeff[s + d + 1] + ((y + 1) ^ d / d - 1) * x$coeff[s + d + 1] * coeff[1] / base::factorial(d)

                coeff[s + d + 1] <- base::factorial(d) * coeff[s + d + 1] / x$coeff[1]
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
        for (k in seq(x$ndirs))
            for (d in seq(x$ord)) {
                s <- (p - 1) * x$ord

                i <- 1
                while (i < d) {
                    s_coeff[s + d + 1] <- s_coeff[s + d + 1] + i * x$coeff[s + i + 1] * c_coeff[s + d - i + 1] / (base::factorial(i) * base::factorial(d - i))
                    c_coeff[s + d + 1] <- c_coeff[s + d + 1] + i * x$coeff[s + i + 1] * s_coeff[s + d - i + 1] / (base::factorial(i) * base::factorial(d - i))
                    i <- i + 1
                }
                s_coeff[s + d + 1] <- s_coeff[s + d + 1] + d * x$coeff[s + d + 1] * c_coeff[1] / base::factorial(d)
                c_coeff[s + d + 1] <- c_coeff[s + d + 1] + d * x$coeff[s + d + 1] * s_coeff[1] / base::factorial(d)

                s_coeff[s + d + 1] <- base::factorial(d - 1) * s_coeff[s + d + 1]
                c_coeff[s + d + 1] <- base::factorial(d - 1) * c_coeff[s + d + 1] * -1
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
                    s_coeff[s + d + 1] <- s_coeff[s + d + 1] + i * x$coeff[s + i + 1] * c_coeff[s + d - i + 1] / (base::factorial(i) * base::factorial(d - i))
                    c_coeff[s + d + 1] <- c_coeff[s + d + 1] + i * x$coeff[s + i + 1] * s_coeff[s + d - i + 1] / (base::factorial(i) * base::factorial(d - i))
                    i <- i + 1
                }
                s_coeff[s + d + 1] <- s_coeff[s + d + 1] + d * x$coeff[s + i + 1] * c_coeff[1] / base::factorial(d)
                c_coeff[s + d + 1] <- c_coeff[s + d + 1] + d * x$coeff[s + i + 1] * s_coeff[1] / base::factorial(d)

                s_coeff[s + d + 1] <- base::factorial(d - 1) * s_coeff[s + d + 1]
                c_coeff[s + d + 1] <- base::factorial(d - 1) * c_coeff[s + d + 1] * -1
            }
    }
    else { #numeric
        return(cos(x))
    }

    return(radx_from(c_coeff, ndirs=x$ndirs))
}
