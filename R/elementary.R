# UTP Algebra for Elementary Functions:
# exp, log, a ^ x, x ^ a, x ^ x, sin, cos

"exp.radx" <- function(x) {
    if(class(x) == "radx") {
        coeff <- rep(0, x$ord + 1)
        coeff[1] <- exp(x$coeff[1])
        for (k in seq(x$ord)) {
            i <- 1
            while (i <= k) {
                coeff[k + 1] <- coeff[k + 1] + i * x$coeff[i + 1] * coeff[k - i + 1] / (base::factorial(i) * base::factorial(k - i))
                i <- i + 1
            }
            coeff[k + 1] <- coeff[k + 1] * base::factorial(k - 1)
        }
    }
    else { #numeric
        return(exp(x))
    }

    return(radx_from(coeff))
}

# natural logarithm (base e)
"log.radx" <- function(x) {
    if(class(x) == "radx") {
        coeff <- rep(0, x$ord + 1)
        coeff[1] <- log(x$coeff[1])
        for (k in seq(x$ord)) {
            i <- 1
            while (i < k) {
                coeff[k + 1] <- coeff[k + 1] + i * coeff[i + 1] * x$coeff[k - i + 1] / (base::factorial(i) * base::factorial(k - i))
                i <- i + 1
            }
            coeff[k + 1] <- (x$coeff[k + 1] - base::factorial(k - 1) * coeff[k + 1]) / x$coeff[1]
        }
    }
    else { #numeric
        return(log(x))
    }

    return(radx_from(coeff))
}

# a ^ x and x ^ a (TODO: x ^ x)
"^.radx" <- function(x, y) {
    if(class(x) == "radx" && class(y) == "radx") { # x ^ x
        return(exp(y * log(x)))
    }
    else if(class(x) == "radx") {
        coeff <- rep(0, x$ord + 1)
        coeff[1] <- x$coeff[1] ^ y
        for (k in seq(x$ord)) {
            i <- 1
            while (i <= k) {
                coeff[k + 1] <- coeff[k + 1] + ((y + 1) ^ i / k - 1) * x$coeff[i + 1] * coeff[k - i + 1] / (base::factorial(i) * base::factorial(k - i))
                i <- i + 1
            }
            coeff[k + 1] <- base::factorial(k) * coeff[k + 1] / x$coeff[1]
        }
    }
    else if(class(y) == "radx") {
        return(exp(y * log(x)))
    }
    else { #numeric
        return(x ^ y)
    }

    return(radx_from(coeff))
}

"sin.radx" <- function(x) {
    if(class(x) == "radx") {
        s_coeff <- rep(0, x$ord + 1)
        c_coeff <- rep(0, x$ord + 1)

        s_coeff[1] <- sin(x$coeff[1])
        c_coeff[1] <- cos(x$coeff[1])
        for (k in seq(x$ord)) {
            i <- 1
            while (i <= k) {
                s_coeff[k + 1] <- s_coeff[k + 1] + i * x$coeff[i + 1] * c_coeff[k - i + 1] / (base::factorial(i) * base::factorial(k - i))
                c_coeff[k + 1] <- c_coeff[k + 1] + i * x$coeff[i + 1] * s_coeff[k - i + 1] / (base::factorial(i) * base::factorial(k - i))
                i <- i + 1
            }
            s_coeff[k + 1] <- base::factorial(k - 1) * s_coeff[k + 1]
            c_coeff[k + 1] <- base::factorial(k - 1) * c_coeff[k + 1] * -1
        }
    }
    else { #numeric
        return(sin(x))
    }

    return(radx_from(s_coeff))
}

"cos.radx" <- function(x) {
    if(class(x) == "radx") {
        s_coeff <- rep(0, x$ord + 1)
        c_coeff <- rep(0, x$ord + 1)

        s_coeff[1] <- sin(x$coeff[1])
        c_coeff[1] <- cos(x$coeff[1])
        for (k in seq(x$ord)) {
            i <- 1
            while (i <= k) {
                s_coeff[k + 1] <- s_coeff[k + 1] + i * x$coeff[i + 1] * c_coeff[k - i + 1] / (base::factorial(i) * base::factorial(k - i))
                c_coeff[k + 1] <- c_coeff[k + 1] + i * x$coeff[i + 1] * s_coeff[k - i + 1] / (base::factorial(i) * base::factorial(k - i))
                i <- i + 1
            }
            s_coeff[k + 1] <- base::factorial(k - 1) * s_coeff[k + 1]
            c_coeff[k + 1] <- base::factorial(k - 1) * c_coeff[k + 1] * -1
        }
    }
    else { #numeric
        return(cos(x))
    }

    return(radx_from(c_coeff))
}
