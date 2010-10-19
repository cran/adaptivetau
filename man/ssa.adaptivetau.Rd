\name{ssa.adaptivetau}
\Rdversion{1.1}
\alias{ssa.adaptivetau}
\alias{adaptivetau}
\title{
Approximate stochastic simulation algorithm
}
\description{
  Implements adaptive tau-leaping approximation for simulating the
  trajectory of a continuous-time Markov process (see reference below).
}
\usage{
ssa.adaptivetau(init.values, transitions, rateFunc, params, tf,
                jacobianFunc = NULL, maxTauFunc = NULL, deterministic = FALSE,
                relratechange = rep(1, length(init.values)),
                tl.params = NA)
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
    length equal to the number of transitions.  Must accept the
    following arguments:
    \itemize{
      \item vector of current values for all state variables
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
  \item{jacobianFunc}{ R-function that returns Jacobian of transition
    rates.  If supplied, enables the use of the implicit tau-leaping
    algorithm (if appropriate; used in stiff systems).  Takes same
    parameters as rateFunc.
  }
  \item{maxTauFunc}{ R-function that returns single real number giving
    the maximum tau leap allowable from the current state.  Only
    executed if the Adaptive tau algorithm wants to leap in a step
    greater than t.params$maxTau (should only need to use if your rate
    functions are dramatically discontinuous).  Takes same parameters as
    rateFunc.
  }
  \item{deterministic}{
    if not FALSE, one-dimensional logical vector of length equal to the
    number of transitions.  Any transitions marked TRUE will be applied
    every timestep using the \bold{expected} degree of change.  Note
    this will almost certainly result in non-integer values for any
    affected state variables.
  }
  \item{relratechange}{
    one-dimensional vector of length equal to the number of state
    variables providing an upperbound to the ratio of amount that any
    transition rate will change givenm a corresponding change in the
    state variable.  In other words, if variable i doubles, can we be
    assured that no transition will \bold{more} than double (ratio of
    1)?  If not, then you need to set this variable to be greater than 1.
  }
  \item{tl.params}{
    Various parameters to the tau-leaping algorithm itself (best
    explained by reading original reference):
    \describe{
      \item{epsilon}{default 0.05; increasing will make bigger leaps
        resulting in potentially more error}
      \item{delta}{default 0.05; how close two symmetric transition
        rates must be before we classify them as in
        partial-equilibrium.  Only applies to the implicit tau routine.
      }
      \item{maxtau}{default Inf; Maximum leap allowed.  Should only need
        to specify if rate functions are drastically non-smooth.
      }
    }
  }
}
\details{
  The initial values, transition matrix & transition rates must all be
  designed such that variables are always non-negative.  The algorithm
  relies on this invariant.

  See reference for details but, in brief the adaptive tau-leaping
  algorithm dynamically switches between three methods for simulating
  events:
  \describe{
    \item{exact}{no approximation -- executes a single transition at a time}
    \item{explicit tau-leaping}{subdivides transitions into those that
      might hit cause a variable to hit 0 (\dQuote{critical}) and those
      that do not.  Duration of time step picked dynamically with the
      goal of maximizing the step while minimized the change in the
      transition rates (the approximation assumes that these rates do
      not change).  Non-critical transitions are advanced by a
      poisson-distributed random variable; critical transitions are
      handled more like the exact algorithm.}
    \item{implicit tau-leaping}{in addition to dividing between critical
      and non-critical, also identifies transitions in quasi-equilbrium
      (reversible pairs of transition that have roughly equal
      forward-backward flow).  Duration of time step picked on basis of
      non-critical, non-equilibrium transitions.  This has the potential
      to greatly increase the timestep size for stiff systems.  Similar
      idea to the explicit method but necesitates solving an implicit
      equation via Newton's method.  Thus you must supply a function to
      calculate the Jacobian to enable this method.}
  }

  All error messages that reference variables or transitions by number
  use 1-based numbering where the first state (or transition) is 1.
}
\value{
  two-dimensional matrix with rows for every timepoint at which the
  rateFunc was evaluated and columns giving the value for every state
  variable at that time.  The first column specifies the time.
}
\references{
  Cao Y, Gillespie DT, Petzold LR, \emph{The Journal of
    Chemical Physics}, 2007
}
\author{
  Philip Johnson
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  For systems with sparse transition matrices, see helper function
  \code{\link{ssa.maketrans}}.  For a pure-R implementation without C++
  underneath, see \code{\link[GillespieSSA]{GillespieSSA}}
}
\examples{
# Lotka-Volterra example
lvrates <- function(x, params, t) {
  with(params, {
    return(c(preygrowth*x[1], ## prey growth rate
             x[1]*x[2]*pred,  ## prey death / predator growth rate
             x[2]*preddeath)) ## predator death rate
  })
}
par=list(preygrowth=10, pred=0.01, preddeath=10);
r=ssa.adaptivetau(c(1000,1000),
                  matrix(c(1,0, -1,1, 0,-1), nrow=2), lvrates, par, tf=1,
                  tl.params=list(epsilon = 0.02))
matplot(r[,"time"], r[,c("x1","x2")], type='l', xlab='Time', ylab='Counts')
legend("topleft", legend=c("prey", "predator"), lty=1:2, col=1:2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
\keyword{ ts }% __ONLY ONE__ keyword per line
