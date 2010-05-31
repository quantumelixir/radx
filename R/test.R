# Test functions for the AD engine
# Grouped into four classes: f(vector or scalar) = (vector or scalar)

#############
# f: R -> R #
#############

t_elem1 <- function (x) {
    sin(x)
}

t_elem2 <- function (x) {
    cos(x)
}

t_elem3 <- function (x) {
    log(x)
}

t_mono <- function (x, n) {
    x ** n
}

# coeff is a vector of coefficients of the polynomial f(x)
t_poly <- function (x, coeff) {
    ret <- 0
    t <- 1
    for (i in coeff) {
        ret <- ret + i * t;
        t <- t * x;
    }
    return(ret)
}

t_quad <- function (x) {
    x*3*x -x*2 + 5
}

t_inv <- function (x) {
    1/x
}

t_simple <- function (x) {
    (x+1)*(x-2)/(x+3)
}

t_alt <- function (x) {
    x - ((4*x + 2)/(x + 3))
}

t_nested <- function (x) {
    exp(sin(exp(cos(x + 2 * x ** 5))))
}

t_example <- function(x) {
    (1 + x + exp(x)) * sin(x)
}

###############
# f: R -> R^n #
###############

###############
# f: R^m -> R #
###############

#################
# f: R^m -> R^n #
#################
