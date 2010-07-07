\name{radx-package}
\alias{radx-package}
\alias{radx}
\docType{package}
\title{
Automatic Differentiation engine in R
}
\description{
Precise numerical evaluation of arbitrary order derivatives for
functions of the form $f:R^m \to R^n$ through both forward
and reverse mode Automatic Differentiation
}
\details{

The main functions are:
\preformatted{
radx        Create objects of class radx. These store the truncated
            Taylor series expansions which are propagated through
            expressions to compute derivatives.

getDeriv    Obtain the derivative from the taylor series coefficients
            of a radx class expression.

radxeval    Driver for computing arbitrary order derivatives of functions.
            Yields a matrix with derivatives of each output function along
            one column ordered by the multi-index corresponding to the
            derivative they represent.

convert2pos Convert a given set of multi-indices (a matrix with one
            multi-index per row) to a vector of positional indices indicating
            the location of the corresponding derivative in the output of radxeval.
            Example: The mixed derivative (d^3 / dx^2 dy^1 dz^0)(f) has the
            multi-index (2, 1, 0) whose positional index is 2. Thus, this
            mixed derivative can be found along the second row of the matrix
            in the output of radxeval.
}
}
\author{
Maintainer: Chidambaram Annamalai <quantumelixir@gmail.com>
}
\references{
    Griewank, A. and Walther, A. (2008) \emph{Evaluating derivatives: principles and techniques of algorithmic differentiation}. Society for Industrial and Applied Mathematics (SIAM)

  Griewank, A. and Utke, J. and Walther, A., (2000) \emph{Evaluating higher derivative tensors by forward propagation of univariate Taylor series}. Journal of Mathematics of Computation, American Mathematical Society.

}
\keyword{ package }
\seealso{
    \code{\link[<numDeriv>:<numDeriv>-package]{numDeriv}}
}
\examples{
\dontrun{
## run demos
demo("univariate")
demo("hessian")

## open the directory with documents
browseURL(paste(system.file(package="radx"), "/doc", sep=""))
}
}
