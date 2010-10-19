\name{ssa.maketrans}
\Rdversion{1.1}
\alias{ssa.maketrans}
\title{
Make transition matrix
}
\description{
  Helper function to (easily) construct sparse transition matrices for
  ssa.adapativetau.
}
\usage{
ssa.maketrans(numVariables, ...)
}
\arguments{
  \item{numVariables}{
    Number of variables in your system.  This will be the number of rows
    in your final matrix.
  }
  \item{...}{
    At least one additional argument is *required*.  Each should be a
    matrix with an arbitrary number of column (each represents a
    separate transition) and an even number of rows.  Going down a
    particular column, entries alternate between specifying the index of
    a variable that will change and the amount by which it will change.
    See examples below.
  }
}
\details{
  Making large transition matrices can be a real pain, particularly if
  the entries in the end will be quite sparse.  The concept of this
  function is simple but easiest way to understand it is probably to
  check out the examples.
}
\value{
  two-dimensional matrix with 'numVariables' rows and columns for each
  transition.
}
\author{
  Philip Johnson
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  You probably only want to use this in preparation for calling
  \code{\link{ssa.adaptivetau}}.
}
\examples{

# trivial Lotka-Volterra example from ssa.adaptivetau
nu = ssa.maketrans(2, rbind(1, +1),
                      rbind(1, -1, 2, +1),
                      rbind(2, -1))

# slightly more complicated SIR model with two distinct susceptible and
# infected variables (think male and female), birth (into S) and death
# (exclusively from I andR)
S = 1:2  # indices of the variables containing number susceptible
I = 3:4  # indices of the variables containing number infected
R = 5    # index of the variable containing number recovered (and immune)
nu = ssa.maketrans(5, rbind(S, +1),
                      rbind(S, -1, I, +1),
                      rbind(I, -1),
                      rbind(I, -1, R, +1),
                      rbind(R, -1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }

