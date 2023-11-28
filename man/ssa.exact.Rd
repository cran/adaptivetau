\name{ssa.exact}
\Rdversion{1.1}
\alias{ssa.exact}
\title{
Exact stochastic simulation algorithm
}
\description{
  Simulates the trajectory of a continuous-time Markov process (aka. the
  Gillespie simulation algorithm).
}
\usage{
ssa.exact(init.values, transitions, rateFunc, params, tf,
          reportTransitions = FALSE)
}
\arguments{
  \item{init.values}{
    Vector of initial values for all variables.  This should be a simple
    one-dimensional vector of real numbers.  \emph{IMPORTANT:} variables
    must never take negative values.
  }
  \item{transitions}{
    Two-dimensional matrix of integers specifying how each state
    variable (rows) should be changed for a given transition (columns).
    Generally this will be a sparse matrix of primarily 1s and -1s.
  }
  \item{rateFunc}{
    R-function that returns instantaneous transition rates for each
    transition in the form a real-valued one-dimensional vector with
    length equal to the number of transitions.  The order of these
    rates must match the order in which transitions were specified in
    the \code{transitions} parameter above.  This function must accept
    the following arguments:
    \itemize{
      \item vector of current values for all state variables (in order
      used in the \code{init.values} argument above)
      \item parameters as supplied in argument to ssa.adaptivetau
      \item single real number giving the current time (all simulations
      start at t=0)
    }
  }
  \item{params}{
    any R variable to be passed as-is to rateFunc, presumably specifying
    useful parameters.
  }
  \item{tf}{
    final time of simulation.  Due to the approximate nature of the
    tau-leaping algorithm, the simulation may overshoot this time somewhat.
  }
  \item{reportTransitions}{default FALSE; whether to include a
    matrix of executed transitions per time interval in the output
    (alongside with the states at each time point).
  }
}
\details{
  This function is supplied as a bonus with the adaptivetau package,
  since the C++ function that underlies this (exact) stochastic
  simulation algorithm is used in the (approximate) adaptive tau
  stochastic simulation as well.
  
  The initial values, transition matrix, and transition rates must all
  be designed such that variables are always non-negative.  The
  algorithm relies on this invariant.

  Consider calling \code{\link[compiler]{enableJIT}(1)} before running
  ssa.exact.  In most circumstances, this should yield some speedup by
  byte-code compiling the rate function.
}
\value{
  Two-dimensional matrix with rows for every timepoint at which the
  rateFunc was evaluated and columns giving the value for every state
  variable at that time.  The first column specifies the time.

  If the \sQuote{reportTransitions} option is used, then a list is
  returned with two elements.  The first element is the
  \sQuote{dynamics} matrix described above.  The second is a
  two-dimensional matrix called \sQuote{transitions} with a row for each
  timepoint and columns giving the number of times each transition was
  executed between the current time and the previous time.
}
\author{
  Philip Johnson
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  This function is a bonus the comes along with the approximate (but
  faster) \code{\link{ssa.adaptivetau}}.  For systems with sparse
  transition matrices, see helper function \code{\link{ssa.maketrans}}.
}
\examples{
## Lotka-Volterra example
lvrates <- function(x, params, t) {
  with(params, {
    return(c(preygrowth*x["prey"],      ## prey growth rate
             x["prey"]*x["pred"]*eat,   ## prey death / predator growth rate
             x["pred"]*preddeath))      ## predator death rate
  })
}
params=list(preygrowth=10, eat=0.01, preddeath=10);
r=ssa.exact(c(prey = 1000, pred = 500),
            matrix(c(1,0, -2,1, 0,-1), nrow=2), lvrates, params, tf=2)
matplot(r[,"time"], r[,c("prey","pred")], type='l', xlab='Time', ylab='Counts')
legend("topleft", legend=c("prey", "predator"), lty=1:2, col=1:2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }

