ssa.adaptivetau <-
function(init.values, transitions, rateFunc, params, tf,
                            jacobianFunc = NULL, maxTauFunc = NULL,
                            deterministic = FALSE,
                            relratechange=rep(1, length(init.values)),
                            tl.params = NA) {
  storage.mode(transitions) = "integer";
  return(.Call('ssa.adaptivetau', PACKAGE='adaptivetau',
               init.values, transitions,
               rateFunc, jacobianFunc, params, tf, deterministic,
               relratechange, tl.params, maxTauFunc))
}

ssa.maketrans <- function(numVariables, ...) {
  trans = list(...)
  if (length(trans) == 0) {
    stop("no transitions passed into ssa.maketrans!")
  }
  ## remember, per convention in the tau-leaping papers, columns
  ## represent transitions and rows represent variables
  
  ## initialize matrix to size (numVariables x num transitions) by
  ## summing over all elements of trans (which may themselves be
  ## matrices)
  m = matrix(0, nrow=numVariables, ncol=sum(sapply(trans,ncol)))

  ## expand out the sparse transitions
  trI = 0
  for (i in 1:length(trans)) {
    x = trans[[i]]
    for (j in 1:ncol(x)) {
      trI = trI + 1
      if (any(x[(1:(nrow(x) %/% 2))*2-1,j] < 1  |
              x[(1:(nrow(x) %/% 2))*2-1,j] > numVariables)) {
        stop("variable index outside valid range (1:numVariables)")
      }
      m[x[(1:(nrow(x) %/% 2))*2-1,j], trI] = x[(1:(nrow(x) %/% 2))*2,j];
    }
  }

  return(m)
}
