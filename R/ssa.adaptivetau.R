ssa.compile <- function(f) {
  if (is.null(f)) return (NULL)
  ## compile function, if R compiler available
  origWarn = options()$warn
  options(warn=-1); #disable warning message if doesn't exist
  haveCompiler = require(compiler, quietly=TRUE)
  options(warn=origWarn)
  if (haveCompiler) {
    compiler::cmpfun(f, options=list(suppressUndefined=TRUE))
  } else {
    f
  }
}

ssa.adaptivetau <-
function(init.values, transitions, rateFunc, params, tf,
         jacobianFunc = NULL, maxTauFunc = NULL,
         deterministic = FALSE,
         relratechange=rep(1, length(init.values)),
         tl.params = NULL) {
  storage.mode(transitions) = "integer";
  return(.Call('simAdaptiveTau', PACKAGE='adaptivetau',
               init.values, transitions,
               ssa.compile(rateFunc), ssa.compile(jacobianFunc),
               params, tf, deterministic,
               relratechange, tl.params, ssa.compile(maxTauFunc)))
}

ssa.exact <-
function(init.values, transitions, rateFunc, params, tf) {
  storage.mode(transitions) = "integer";
  return(.Call('simExact', PACKAGE='adaptivetau',
               init.values, transitions, ssa.compile(rateFunc), params, tf))
}

ssa.maketrans <- function(variables, ...) {
  trans = list(...)
  if (length(trans) == 0) {
    stop("no transitions passed into ssa.maketrans!")
  }
  if (length(variables) == 1  &&  is.numeric(variables)) {
    numVariables = variables;
  } else if (is.character(variables)) {
    numVariables = length(variables)
  } else {
    stop("Cannot deduce number of variables -- ssa.maketrans requires ",
         "either a vector of variable names or the number of variables")
  }

  ## initialize matrix to size (numVariables x num transitions) by
  ## summing over all elements of trans (which may themselves be
  ## matrices)
  ## per convention in the tau-leaping papers, columns represent
  ## transitions and rows represent variables
  m = matrix(0, nrow=numVariables, ncol=sum(sapply(trans,ncol)))
  if (is.character(variables)) {
    rownames(m) = variables;
  }

  ## expand out the sparse transitions
  trI = 0
  for (i in 1:length(trans)) {
    x = trans[[i]]
    idx = (1:(nrow(x) %/% 2))*2-1 #indices of variables
    mag = idx + 1                 #indices of magnitudes
    for (j in 1:ncol(x)) {
      trI = trI + 1
      if (is.numeric(x)) {
        if (any(x[idx,j] < 1  |  x[idx,j] > numVariables)) {
          stop("variable index outside valid range (1:numVariables)")
        }
        m[x[idx,j], trI] = x[mag,j];
      } else if (is.character(x)) {
        if (any(!(x[idx,j] %in% rownames(m)))) {
          stop("unknown variable(s): ",
               paste(x[idx,j][!(x[idx,j] %in% rownames(m))], collapse=", "))
        }
        m[x[idx,j], trI] = as.double(x[mag,j])
      } else {
        stop("transitions passed to ssa.maketrans must be integer or character")
      }
    }
  }

  return(m)
}
