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
ssa.maketrans(variables, ...)
}
\arguments{
  \item{variables}{
    Either the number of variables in your system OR a vector of strings
    specifying the names of the variables (in which case the length of
    this vector is the number of variables).  The number of variables
    equals the number of rows in the returned transition matrix.
  }
  \item{...}{
    At least one additional argument is *required*.  Each should be a
    matrix with an arbitrary number of columns (each represents a
    separate transition) and an even number of rows.  Going down a
    particular column, entries alternate between specifying the variable
    that will change and the amount by which it will change.  Variables
    be specified either by numerical index (1,2,...,numVariables) or
    character name (if names supplied in the parameter above).  See
    examples below.
  }
}
\details{
  Making large transition matrices can be a real pain, particularly if
  the entries in the end will be quite sparse.  The concept of this
  function is simple but the easiest way to understand it is probably to
  check out the examples.
}
\value{
  two-dimensional matrix with rows for each variable and columns for each
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

## trivial Lotka-Volterra example from ssa.adaptivetau
nu = ssa.maketrans(2, #number of variables
                   rbind(1, +1),
                   rbind(1, -1, 2, +1),
                   rbind(2, -1))

## slightly more complicated SIR epidemiological model with two distinct
## susceptible and infected variables (think male and female), birth
## (into S) and death (exclusively from I and R)
nu = ssa.maketrans(c("Sm", "Sf", "Im", "If", "R"), #list of variable names
                   rbind(c("Sm","Sf"), +1),
                   rbind(c("Sm","Sf"), -1, c("Im","If"), +1),
                   rbind(c("Im","If"), -1),
                   rbind(c("Im","If"), -1, "R", +1),
                   rbind("R", -1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }

