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
                jacobianFunc = NULL, maxTauFunc = NULL,
                deterministic = NULL, halting = NULL,
                relratechange = rep(1, length(init.values)),
                tl.params = NULL)
}
\arguments{
  \item{init.values}{
    Vector of initial values for all variables.  This should be a simple
    one-dimensional vector of real numbers.  \emph{IMPORTANT:} variables
    must never take negative values.
  }
  \item{transitions}{
    One of two possible data types:
    \itemize{
      \item A list with length equal to the number of transitions.  Each
      element of the list should be a vector specifying a transition
      (i.e., which state(s) change and by how much).  Each entry in the
      vector needs a name (specifying which state variable to change,
      either by name or index) and a value (specifying the amount by
      which this variable will change).
      \item A two-dimensional matrix of integers specifying how each state
      variable (rows) should be changed for a given transition (columns).
      Generally this will be a sparse matrix of primarily 1s and -1s,
      which can make this structure inefficient.
    }
    See the example below for details as well as
    \code{\link{ssa.maketrans}} or the vignette accompanying this package.
  }
  \item{rateFunc}{
    R function that returns instantaneous transition rates for each
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
  \item{jacobianFunc}{ R function that returns Jacobian of transition
    rates.  If supplied, enables the use of the implicit tau-leaping
    algorithm (if appropriate; used in stiff systems).  In contrast to
    typical definition of Jacobian (but consistent with the transition
    matrix parameter above), columns represent transitions (functions)
    and rows represent state variables (different partial derivative
    variables). jacobianFunc takes the same parameters as rateFunc.
    \bold{NOTE:} the implicit tau-leaping algorithm will reduce the observed
    variance for fast-changing state variables, although the mean will
    be correct.
  }
  \item{maxTauFunc}{ R function that returns single real number giving
    the maximum tau leap allowable from the current state.  Only
    executed if the adaptive tau algorithm wants to leap in a step
    greater than tl.params$maxtau (should only need to use if your rate
    functions have a dramatically discontinuous first derivative).
    Takes same parameters as rateFunc.
  }
  \item{deterministic}{
    Specify transitions to be treated as deterministic.  If not NULL,
    either a TRUE/FALSE vector of length equal to the number of
    transitions -or- a integer vector of transition numbers to be
    treated as deterministic (i.e. applying \code{\link{which}} to the
    former).  Deterministic transitions will be applied every timestep
    using the \bold{expected} degree of change.  Note this will almost
    certainly result in non-integer values for any affected state
    variables.  Also, note that this is still an approximation -- not a
    true numerical ODE solution -- and at least one transition must be
    stochastic for this to work at all.
  }
  \item{halting}{
    Specify transitions which, when executed, will cause the simulation
    to halt.  If not NULL, either a TRUE/FALSE vector of length equal to the
    number of transitions -or- a integer vector of transition numbers to
    be treated as halting (i.e. applying \code{\link{which}} to
    the former).  If halting transitions are specified, then the return
    value is different (since we want to know which transition halted us).
  }
  \item{relratechange}{
    one-dimensional vector of length equal to the number of state
    variables providing an upperbound to the ratio of amount that any
    transition rate will change given a corresponding change in the
    state variable.  In other words, if variable i doubles, can we be
    assured that no transition will \bold{more} than double (ratio of
    1)?  If not, then you need to set this variable to be greater than 1.
  }
  \item{tl.params}{
    List of various parameters to the tau-leaping algorithm itself (best
    explained by reading original reference):
    \describe{
      \item{epsilon}{default 0.05; increasing will make bigger leaps
        resulting in potentially more error}
      \item{delta}{default 0.05; how close two symmetric transition
        rates must be before we classify them as in
        partial-equilibrium.  Only applies to the implicit tau routine.
      }
      \item{maxtau}{default Inf; maximum leap allowed.  Should only need
        to specify if rate functions are drastically non-smooth.
      }
      \item{extraChecks}{default TRUE; whether to perform a battery of
        validity checks on the state variables and rates after every
        call to "rateFunc."  Slows down implementation, but unless
        you are SURE you have designed your transition matrix and rate
        function correctly, probably best to keep these checks.
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
  use 1-based numbering where the first variable (or transition) is 1.

  Consider calling \code{\link[compiler]{enableJIT}(1)} before running
  ssa.adaptivetau.  In most circumstances, this should yield some
  speedup by byte-code compiling the rate function.
}
\value{
  If no halting transitions are specified, then a two-dimensional matrix
  with rows for every timepoint at which the rateFunc was evaluated and
  columns giving the value for every state variable at that time.  The
  first column specifies the time.

  If halting transitions are specified, then a list with two elements.
  The first (\sQuote{dynamics}) is the same two-dimensional matrix as
  without halting transitions.  The second (\sQuote{haltingTransition})
  gives the number of the transition that halted the simulation, or NA
  otherwise.
}
\note{
  Development of this package was supported in part by National Science
  Foundation award DBI-0906041 and National Institute of Health award
  K99-GM104158.  Any opinions, findings, and conclusions or
  recommendations expressed in this material are those of the author and
  do not necessarily reflect the views of the NSF or NIH.
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
  underneath, see \code{\link[GillespieSSA]{GillespieSSA}}.  Also,
  \code{\link{ssa.exact}} exposes a R interface to the C++
  implementation of the exact, non-approximate simulation algorithm
  (sometimes known as "Gillespie").
}
\examples{

## Simple Lotka-Volterra example

# We have three potential transitions:
transitions = list(c(prey = +1), # prey grow
                   c(prey = -2, pred = +1), # predation
                   c(pred = -1)) # predator dies

# Function to calculate transition rates, given variables and parameters
lvrates <- function(x, params, t) {
  return(c(params$preygrowth*x["prey"],     # rate of prey growing
           x["prey"]*x["pred"]*params$eat,  # rate of predation
           x["pred"]*params$preddeath))     # rate of predators dying
}

# Set the Lotka-Volterra parameters
params = list(preygrowth=10, eat=0.01, preddeath=10);

# Set the random seed (only necessary if you want to reproduce results)
set.seed(4)

# Perform the stochastic simulation!
r=ssa.adaptivetau(c(prey = 1000, pred = 500),
                  transitions, lvrates, params, tf=12)

# Plot the results
matplot(r[,"time"], r[,c("prey","pred")], type='l', xlab='Time',
        ylab='Counts (log scale)', log='y')
legend("bottomleft", legend=c("prey", "predator"), lty=1:2, col=1:2)


# However, if you are interested in very fine-scale variance, perhaps the
# default parameters smooth too much.  Try reducing the tl.param epsilon
# to see a better approximation:
s=ssa.adaptivetau(c(prey = 1000, pred = 500),
                  transitions, lvrates, params, tf=12,
                  tl.params = list(epsilon=0.01)) # reduce "epsilon"
par(mfrow=c(2,1));
matplot(r[r[,"time"]<2,"time"], r[r[,"time"]<2,c("prey","pred")],
        type='l', xlab='Time', ylab='Counts', main="Original (epsilon=default)")
matplot(s[s[,"time"]<2,"time"], s[s[,"time"]<2,c("prey","pred")],
        type='l', xlab='Time', ylab='Counts', main="Fine-scale (epsilon=0.01)")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
\keyword{ ts }% __ONLY ONE__ keyword per line
