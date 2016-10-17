\name{NNLS}
\alias{NNLS}
\title{
Sequential Coordinate-wise Algorithm for the Non-negative Least Squares Problem
}
\description{
Consider the linear system \eqn{\bold{A} x = b} where \eqn{\bold{A} \in  R\textsuperscript{m x n}}, \eqn{x \in R\textsuperscript{n}}, and \eqn{b \in R\textsuperscript{m}}.  The technique of least squares proposes to compute \eqn{x} so that the sum of squared residuals is minimized.  \code{NNLS} solves the least squares problem \eqn{\min{||\bold{A} x = b||\textsuperscript{2}}} subject to the constraint \eqn{x \ge 0}.  This implementation of the Sequential Coordinate-wise Algorithm uses a sparse input matrix \eqn{\bold{A}}, which makes it efficient for large sparse problems.
}
\usage{
NNLS(A,
     b,
     precision = sqrt(.Machine$double.eps),
     processors = 1,
     verbose = TRUE)
}
\arguments{
  \item{A}{
List representing the sparse matrix with integer components i and j, numeric component x.  The fourth component, dimnames, is a list of two components that contains the names for every row (component 1) and column (component 2).
}
  \item{b}{
Numeric matrix for the set of observed values.  (See details section below.)
}
  \item{precision}{
The desired accuracy.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
The input \eqn{b} can be either a matrix or a vector of numerics.  If it is a matrix then it is assumed that each column contains a set of observations, and the output \eqn{x} will have the same number of columns.  This allows multiple NNLS problems using the same \eqn{\bold{A}} matrix to be solved simultaneously, and greatly accelerates computation relative to solving each sequentially.
}
\value{
A list of two components:
  \item{x}{
The matrix of non-negative values that best explains the observed values given by \eqn{b}.
}
  \item{res}{
A matrix of residuals given by \eqn{\bold{A} x - b}.
}
}
\references{
Franc, V., et al. (2005). Sequential coordinate-wise algorithm for the non-negative least squares problem.  Computer Analysis of Images and Patterns, 407-414.
}
\seealso{
\code{\link{Array2Matrix}}, \code{\link{DesignArray}}
}
\examples{
# unconstrained least squares:
A <- matrix(c(1, -3, 2, -3, 10, -5, 2, -5, 6), ncol=3)
b <- matrix(c(27, -78, 64), ncol=1)
x <- solve(crossprod(A), crossprod(A, b))

# Non-negative least squares:
w <- which(A > 0, arr.ind=TRUE)
A <- list(i=w[,"row"], j=w[,"col"], x=A[w],
          dimnames=list(1:dim(A)[1], 1:dim(A)[2]))
x_nonneg <- NNLS(A, b)

# compare the unconstrained and constrained solutions:
cbind(x, x_nonneg$x)

# the input value "b" can also be a matrix:
b2 <- matrix(b, nrow=length(b), ncol=2) # repeat b in two columns
x_nonneg <- NNLS(A, b2) # solution is repeated in two output columns
}
