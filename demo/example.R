# define a two input, two ouput vector function
f <- function(x, y) {
c(1 - 2*y + sin(4*pi*y) - x, y - sin(2*pi*x)/2)
}

# find its derivatives upto the first order about the point (x=1.1,y=2.2)
radxeval(f, c(1.1,2.2), 1)

# find its derivatives upto the second order about the point (x=1.1,y=2.2)
radxeval(f, c(1.1,2.2), 2)

# The only problem is reading off the results as the output is in a compressed tensor derivative notation. 
# We use convert2pos to convert a multi-index corresponding to a partial derivative to an index in the array output by radxeval
# Basically the (i,j) element in the ouput matrix stores the ith derivative of the jth output function. By ith derivative I mean the partial derivative corresponding to the ith multi-index. To simplify these issues, we can now: 
hess <- radxeval(f, c(1.1,2.2), 2)

# extract the d^2/dx^0dy^2 derivative (which is where the 0, 2 come from) of the first output function (which is where the 1 comes from)
hess[convert2pos(c(0,2)), 1]
