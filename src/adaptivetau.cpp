/* $Id$
    --------------------------------------------------------------------------
    C++ implementation of the "adaptive tau-leaping" algorithm described by
    Cao Y, Gillespie DT, Petzold LR. The Journal of Chemical Physics (2007).
    Author: Philip Johnson <plfjohnson@emory.edu>


    Copyright (C) 2010 Philip Johnson

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


    If building library outside of R package (i.e. for debugging):
        R CMD SHLIB adaptivetau.cpp
    --------------------------------------------------------------------------
*/

#include <iostream>
#include <vector>
#include <limits>
#include <sstream>
#include <stdexcept>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>

#include "Rwrappers.h"

using namespace std;

enum EStepType {
    eExact = 0,
    eExplicit,
    eImplicit
};

const bool debug = false;

class CEarlyExit : public runtime_error {
public:
    CEarlyExit(const string &w) : runtime_error(w) {}
};

//use below rather than R's "error" directly (which will not free memory, etc.)
#ifdef throwError
#undef throwError
#endif
#define throwError(e) { ostringstream s; s << e; throw runtime_error(s.str()); }
#ifdef throwEarlyExit
#undef throwEarlyExit
#endif
#define throwEarlyExit(e) { ostringstream s; s << e << "; results returned only up until this point"; throw CEarlyExit(s.str());}

// Functions below are a hack suggested by Simon Urbanek (although
// he "would not recommend for general use") to check if the user has
// asked to interrupt execution.  The problem with calling
// R_CheckUserInterrupt directly is it longjmps and we don't have
// a chance to free memory from the heap.
// http://tolstoy.newcastle.edu.au/R/e13/devel/11/04/1049.html
// See extern "C" definition of chkIntFn above.
extern "C" void chkIntFn(void*) { R_CheckUserInterrupt(); }
bool AdaptiveTauCheckUserInterrupt(void) {
    return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

class CStochasticEqns {
public:
    CStochasticEqns(SEXP initVal, SEXP nu,
                    SEXP rateFunc, SEXP rateJacobianFunc,
                    SEXP params, double* changeBound, SEXP maxTauFunc,
                    SEXP detTrans, SEXP haltTrans, SEXP reportTransitions) {
        // copy initial values into new vector (keeping in SEXP vector
        // allows easy calling of R function to calculate rates)
        m_NumStates = Rf_length(initVal);
        SEXP x;
        x = PROTECT(Rf_allocVector(REALSXP, m_NumStates));//protected until ~CStochasticEqns
        Rf_copyVector(x, Rf_coerceVector(initVal,REALSXP));
        if (Rf_isNull(Rf_getAttrib(initVal, R_NamesSymbol))) {
            m_VarNames = NULL;
        } else {
            SEXP namesO = PROTECT(Rf_getAttrib(initVal, R_NamesSymbol));
            m_VarNames = PROTECT(Rf_allocVector(STRSXP, Rf_length(namesO)));
            Rf_copyVector(m_VarNames, namesO);
            Rf_setAttrib(x, R_NamesSymbol, m_VarNames);
            UNPROTECT(2);
            PROTECT(m_VarNames);
        }
        m_X = REAL(x);

        // copy Nu matrix into my own sparse matrix data structure
        if (Rf_isMatrix(nu)) { //old matrix data structure
            CRMatrix<int> mat(PROTECT(Rf_coerceVector(nu,INTSXP)));
            m_Nu.resize(mat.ncol());
            for (int i = 0;  i < mat.nrow();  ++i) {
                for (int j = 0;  j < mat.ncol();  ++j) {
                    if (mat(i,j) != 0) {
                        SChange s;
                        s.m_State = i; s.m_Mag = mat(i,j);
                        m_Nu[j].push_back(s);
                    }
                }
            }
            UNPROTECT(1);
        } else { //list (newer, sparse data structure)
            CRList list(nu);
            m_Nu.resize(list.size());
            for (unsigned int j = 0;  j < list.size();  ++j) {
                if (!Rf_isInteger(list[j])  &&  !Rf_isReal(list[j])) {
                    throwError("the sparse transition matrix representation "
                               "must be a list of either integer or double "
                               "vectors.");
                }
                const CRVector<int> trans(PROTECT(Rf_coerceVector(list[j],INTSXP)));
                UNPROTECT(1);
                m_Nu[j].resize(trans.size());
                for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                    short int state = -1;
                    const char *stateStr = trans.GetName(i);
                    if (strcmp(stateStr, "") == 0) {
                        throwError("transition matrix contains values without "
                                   "a corresponding state variable.");
                    }
                    if (m_VarNames != NULL) {
                        for (state = 0;  state < Rf_length(m_VarNames)  &&
                                 strcmp(CHAR(STRING_PTR(m_VarNames)[state]),
                                        stateStr) != 0;
                             ++state);
                    }
                    if (state < 0  ||  state >= (int) m_NumStates) {
                        istringstream iss(stateStr);
                        iss >> state;
                        if (!iss  ||  !iss.eof()) {
                            state = -1;
                        } else {
                            --state; //switch from 1-based to 0-based
                        }
                    }
                    if (state < 0  ||  state >= (int) m_NumStates) {
                        throwError("transition matrix references non-existent "
                                   "state variable '" << stateStr << "'");
                    }
                    m_Nu[j][i].m_State = state;
                    m_Nu[j][i].m_Mag = trans[i];
                }
            }
        }
        m_TransCats.resize(m_Nu.size(), eNormal);

        // potentially flag some transitions as "deterministic"
        x_SetCat(detTrans, eDeterministic);
        if (m_TransByCat[eDeterministic].size() == m_TransCats.size()) {
            throwError("At least one transition must be stochastic (all "
                       "transitions are currently flagged as "
                       "deterministic).");
        }

        // potentially flag some transitions as "halting" (which are
        // always critical)
        x_SetCat(haltTrans, eHalting);
        m_TransByCat[eCritical] = m_TransByCat[eHalting];

        // needed for ITL
        x_IdentifyBalancedPairs();
        x_IdentifyRealValuedVariables();

        // prepare R function for evaluation by setting up arguments
        // (current X values, parameters, current time)
        SEXP s_time;
        s_time = PROTECT(Rf_allocVector(REALSXP, 1));//protected until ~CStochasticEqns
        m_T = REAL(s_time);
        m_RateFunc = PROTECT(Rf_lang4(rateFunc, x, params, s_time));//protected until ~CStochasticEqns
        if (!rateJacobianFunc  ||  Rf_isNull(rateJacobianFunc)) {
            m_RateJacobianFunc = NULL;
        } else {
            PROTECT(m_RateJacobianFunc = Rf_lang4(rateJacobianFunc, x,
                                               params, s_time)); //protected until ~CStochasticEqns
        }
        m_Rates = NULL;

        //default parameters to adaptive tau leaping algorithm
        m_Epsilon = 0.05;
        m_Ncritical = 10;
        m_Nstiff = 100;
        m_ExactThreshold = 10;
        m_Delta = 0.05;
        m_NumExactSteps[eExact] = 100;
        m_NumExactSteps[eExplicit] = 100;
        m_NumExactSteps[eImplicit] = 10;
        m_ITLConvergenceTol = 0.01;
        m_MaxTau = numeric_limits<double>::infinity();
        m_MaxSteps = 0; // special case 0 == no limit

        //useful additional parameters
        m_RecordTransitionTimeSeries = LOGICAL(reportTransitions)[0];
        m_ExtraChecks = true;
        m_VerboseTracing = 0;
        m_RateChangeBound = changeBound;
        if (!maxTauFunc  ||  Rf_isNull(maxTauFunc)) {
            m_MaxTauFunc = NULL;
        } else {
            PROTECT(m_MaxTauFunc = Rf_lang4(maxTauFunc, x, params, s_time));//protected until ~CStochasticEqns
        }

        //check initial conditions to make sure legit
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            if (m_X[i] < 0) {
                throwError("initial value for variable " << i+1 <<
                           " must be positive (currently " << m_X[i] << ")");
            }
            if (!m_RealValuedVariables[i]  &&
                (m_X[i] - Rf_ftrunc(m_X[i]) > 1e-5)) {
                if (m_VarNames != NULL) {
                    throwError("initial value for variable " << i+1 <<
                               " ('"<<CHAR(STRING_PTR(m_VarNames)[i])<<"') " <<
                               "must be an integer (currently " << m_X[i] << ")");
                } else {
                    throwError("initial value for variable " << i+1 <<
                               " must be an integer (currently " << m_X[i] << ")");
                }
            }
        }
        
        *m_T = 0;
        m_LastTransition = -1;
        m_PrevStepType = eExact;
        GetRNGstate();
    }
    ~CStochasticEqns(void) {
        int cnt = 4;
        if (m_RateJacobianFunc != NULL) {
            ++cnt;
        }
        if (m_Rates != NULL) {
            ++cnt;
        }
        if (m_MaxTauFunc != NULL) {
            ++cnt;
        }
        UNPROTECT(cnt);
    }
    void SetTLParams(SEXP list) {
        SEXP names = PROTECT(Rf_getAttrib(list, R_NamesSymbol));
        try {
            for (int i = 0;  i < Rf_length(names);  ++i) {
                if (strcmp("epsilon", CHAR(STRING_PTR(names)[i])) == 0) {
                    if (!Rf_isReal(VECTOR_ELT(list, i))  ||
                        Rf_length(VECTOR_ELT(list, i)) != 1) {
                        throwError("invalid value for parameter '" <<
                                   CHAR(STRING_PTR(names)[i]) << "'");
                    }
                    m_Epsilon = REAL(VECTOR_ELT(list, i))[0];
                } else if (strcmp("delta", CHAR(STRING_PTR(names)[i])) == 0) {
                    if (!Rf_isReal(VECTOR_ELT(list, i))  ||
                        Rf_length(VECTOR_ELT(list, i)) != 1) {
                        throwError("invalid value for parameter '" <<
                                   CHAR(STRING_PTR(names)[i]) << "'");
                    }
                    m_Delta = REAL(VECTOR_ELT(list, i))[0];
                } else if (strcmp("maxtau", CHAR(STRING_PTR(names)[i])) == 0) {
                    if (!Rf_isReal(VECTOR_ELT(list, i))  ||
                        Rf_length(VECTOR_ELT(list, i)) != 1) {
                        throwError("invalid value for parameter '" <<
                                   CHAR(STRING_PTR(names)[i]) << "'");
                    }
                    m_MaxTau = REAL(VECTOR_ELT(list, i))[0];
                } else if (strcmp("extraChecks",
                                  CHAR(STRING_PTR(names)[i])) == 0) {
                    if (!Rf_isLogical(VECTOR_ELT(list, i))  ||
                        Rf_length(VECTOR_ELT(list, i)) != 1) {
                        throwError("invalid value for parameter '" <<
                                   CHAR(STRING_PTR(names)[i]) << "'");
                    }
                    m_ExtraChecks = LOGICAL(VECTOR_ELT(list, i))[0];
                } else if (strcmp("verbose",
                                  CHAR(STRING_PTR(names)[i])) == 0) {
                    if (!Rf_isInteger(VECTOR_ELT(list, i))  ||
                        Rf_length(VECTOR_ELT(list, i)) != 1) {
                        throwError("invalid value for parameter '" <<
                                   CHAR(STRING_PTR(names)[i]) << "'");
                    }
                    m_VerboseTracing = INTEGER(VECTOR_ELT(list, i))[0];
                } else {
                    Rf_warning("ignoring unknown parameter '%s'",
                            CHAR(STRING_PTR(names)[i]));
                }
            }
        } catch (...) {
            UNPROTECT(1);
            throw;
        }
        UNPROTECT(1);
    }

    void EvaluateATLUntil(double tF) {
        unsigned int c = 0;
        //add initial conditions to time series
        m_TimeSeries.push_back(STimePoint(0, m_X, m_NumStates));
        if (m_RecordTransitionTimeSeries) {
            m_ExecutedTransitions.assign(m_Nu.size(), 0);
            m_TransitionTimeSeries.push_back(m_ExecutedTransitions);
        }
        //main loop
        while (*m_T < tF  &&  (m_MaxSteps == 0 || c < m_MaxSteps)  &&
               (m_LastTransition < 0  ||
                m_TransCats[m_LastTransition] != eHalting)) {
            x_UpdateRates();
            x_SingleStepATL(tF);
            if (++c % 10 == 0  &&  AdaptiveTauCheckUserInterrupt()) {
                throwEarlyExit("simulation interrupted by user at time " << *m_T
                               << " after " << c << " time steps.");
            }
        }
        //save RNG state back to R (could also do in destructor, but
        //no harm in extra calls to PutRNGstate and avoids potential
        //problems with PROTECTing return value from GetTimeSeriesSEXP
        PutRNGstate();
    }
    void EvaluateExactUntil(double tF) {
        unsigned int c = 0;
        //add initial conditions to time series
        m_TimeSeries.push_back(STimePoint(0, m_X, m_NumStates));
        if (m_RecordTransitionTimeSeries) {
            m_ExecutedTransitions.assign(m_Nu.size(), 0);
            m_TransitionTimeSeries.push_back(m_ExecutedTransitions);
        }
        m_LastTransition = -1;
        //main loop
        while (*m_T < tF  &&  (m_MaxSteps == 0 || c < m_MaxSteps)  &&
               (m_LastTransition < 0  ||
                m_TransCats[m_LastTransition] != eHalting)) {
            x_UpdateRates();
            x_SingleStepExact(tF);
            if (++c % 10 == 0  &&  AdaptiveTauCheckUserInterrupt()) {
                throwEarlyExit("simulation interrupted by user at time " << *m_T
                               << " after " << c << " time steps.");
            }
        }
        //save RNG state back to R (could also do in destructor, but
        //no harm in extra calls to PutRNGstate and avoids potential
        //problems with PROTECTing return value from GetTimeSeriesSEXP
        PutRNGstate();
    }

    SEXP GetResult(void) const {
        if (m_TransByCat[eHalting].size() == 0) {
            if (m_RecordTransitionTimeSeries) {
                CRList res(2);
                PROTECT(res);
                res.SetSEXP(0, PROTECT(GetTimeSeriesSEXP()), "dynamics");
                res.SetSEXP(1, PROTECT(GetTransitionTimeSeriesSEXP()), "transitions");
                UNPROTECT(3);
                return res;
            } else {
                return GetTimeSeriesSEXP();
            }
        } else {
            CRList res(m_RecordTransitionTimeSeries ? 3 : 2);
            PROTECT(res);
            res.SetSEXP(0, PROTECT(GetTimeSeriesSEXP()), "dynamics");
            CRVector<int> lastTrans(1);
            lastTrans[0] = (m_LastTransition < 0  ||
                            m_TransCats[m_LastTransition] != eHalting) ?
                NA_INTEGER : m_LastTransition+1;
            res.SetSEXP(1, lastTrans, "haltingTransition");
            if (m_RecordTransitionTimeSeries) {
                res.SetSEXP(2, PROTECT(GetTransitionTimeSeriesSEXP()), "transitions");
                UNPROTECT(1);
            }
            UNPROTECT(2);
            return res;
        }
    }

    SEXP GetTimeSeriesSEXP(void) const {
        SEXP res;
        PROTECT(res = Rf_allocMatrix(REALSXP, m_TimeSeries.size(), m_NumStates+1));
        double *rvals = REAL(res);
        for (unsigned int t = 0;  t < m_TimeSeries.size();  ++t) {
            rvals[t] = m_TimeSeries[t].m_T;
            for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                rvals[(i+1) * m_TimeSeries.size() + t] = m_TimeSeries[t].m_X[i];
            }
        }

        SEXP dimnames, colnames;
        PROTECT(dimnames = Rf_allocVector(VECSXP, 2));
        PROTECT(colnames = Rf_allocVector(VECSXP, m_NumStates+1));
        SET_VECTOR_ELT(dimnames, 1, colnames);
        SET_VECTOR_ELT(colnames, 0, Rf_mkChar("time"));
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            if (m_VarNames  &&  (unsigned int)Rf_length(m_VarNames) > i) {
                SET_VECTOR_ELT(colnames, i+1,
                               STRING_PTR(m_VarNames)[i]);
            } else {
                ostringstream oss;
                oss << "x" << i+1;
                SET_VECTOR_ELT(colnames, i+1, Rf_mkChar(oss.str().c_str()));
            }
        }
        Rf_setAttrib(res, R_DimNamesSymbol, dimnames);

        UNPROTECT(3);
        return res;
    }
    SEXP GetTransitionTimeSeriesSEXP(void) const {
        SEXP res;
        PROTECT(res = Rf_allocMatrix(REALSXP, m_TimeSeries.size(), m_Nu.size()));
        double *rvals = REAL(res);
        for (unsigned int t = 0;  t < m_TimeSeries.size();  ++t) {
            for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
                rvals[j * m_TimeSeries.size() + t] = m_TransitionTimeSeries[t][j];
            }
        }
        UNPROTECT(1);
        return res;
    }

protected:
    enum ETransCat {
        eNormal = 0,
        eCritical,
        eDeterministic,
        eHalting
    };
    typedef vector<ETransCat> TTransCats;
    typedef vector<int> TTransList;
    typedef vector<pair<unsigned int, unsigned int> > TBalancedPairs;
    typedef vector<bool> TBools;
    typedef double* TStates;
    typedef double* TRates;
    struct SChange {
        short int m_State;
        short int m_Mag;
    };
    typedef vector< vector<SChange> > TTransitions;
    struct STimePoint {
        STimePoint(double t, double *x, int n) {
            m_T = t;
            m_X = new double[n];
            memcpy(m_X, x, n*sizeof(double));
        }
        double m_T;
        double *m_X;
    };
    class CTimeSeries : public vector<STimePoint> {
    public:
        ~CTimeSeries(void) {
            for (iterator i = begin();  i != end();  ++i) {
                delete[] i->m_X; i->m_X = NULL;
            }
        }
    };
    typedef vector<double> TExecutedTransitions;
    typedef vector<TExecutedTransitions> TTransitionTimeSeries;

protected:
    void x_IdentifyBalancedPairs(void);
    void x_IdentifyRealValuedVariables(void);
    void x_SetCat(SEXP trans, ETransCat cat);

    void x_AdvanceDeterministic(double deltaT, bool clamp = false);
    void x_SingleStepExact(double tf);
    void x_SingleStepETL(double tau);
    void x_SingleStepITL(double tau);
    void x_SingleStepATL(double tf);

    void x_UpdateRates(void) {
        if (m_ExtraChecks) {
            for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                if (m_X[i] < 0) {
                    throwError("negative variable: " << i+1 << " is " <<
                               m_X[i] << " (check rate function "
                               "and/or transition matrix)");
                } else if (ISNAN(m_X[i])) {
                    throwError("NaN variable: " << i+1 << " is " <<
                               m_X[i] << " (check rate function "
                               "and/or transition matrix)");
                }
            }
        }

        //sigh.  Cover ourselves in case rateFunc code messes with
        //RNGs (which really should NOT be the case -- we rely on the
        //rate functions being deterministic!), but ran into this
        //problem when using a Rcpp rate function (which adds
        //arbitrary calls to Get/Put RNG).
        PutRNGstate(); 

        // make sure our rates are protected!
        if (m_Rates != NULL) { 
            UNPROTECT(1);
            m_Rates = NULL;
        }
        SEXP res = PROTECT(Rf_eval(m_RateFunc, R_EmptyEnv));
        m_Rates = REAL(res);

        if ((unsigned int) Rf_length(res) != m_Nu.size()) {
            throwError("invalid rate function -- returned number of rates ("
                       << Rf_length(res) << ") is not the same as specified by "
                       "the transition matrix (" << m_Nu.size() << ")!");
        }
        if (m_ExtraChecks) {
            for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
                if (ISNAN(m_Rates[j])) {
                    throwError("invalid rate function -- rate for transition "
                               << j+1 << " is not a number (NA/NaN)! (check "
                               "for divison by zero or similar)");
                }
                if (m_Rates[j] < 0) {
                    throwError("invalid rate function -- rate for transition "
                               << j+1 << " is negative!");
                }
            }
        }
    }
    double* x_CalcJacobian(void) {
        SEXP res = Rf_eval(m_RateJacobianFunc, R_EmptyEnv);
        if (!Rf_isMatrix(res)) {
            throwError("invalid Jacobian function -- should return a " <<
                       m_NumStates << " by " << m_Nu.size() << " matrix");
        }
        unsigned int nrow = INTEGER(Rf_getAttrib(res, R_DimSymbol))[0];
        unsigned int ncol = INTEGER(Rf_getAttrib(res, R_DimSymbol))[1];
        if (nrow != m_NumStates  ||  ncol != m_Nu.size()) {
            throwError ("invalid Jacobian function -- returned a " << nrow
                        << " by " << ncol << " matrix instead of the expected "
                        << m_NumStates << " by " << m_Nu.size() <<
                        " (variables by transitions)");
        }
        return REAL(res);
    }
    double x_CalcUserMaxTau(void) {
        if (!m_MaxTauFunc) { throwError("logic error at line " << __LINE__) }
        SEXP res = Rf_eval(m_MaxTauFunc, R_EmptyEnv);
        if (Rf_length(res) != 1  || !Rf_isReal(res)) {
            throwError("invalid return value from maxTau function (should be "
                       "a single real number)");
        }
        return REAL(res)[0];
    }

    unsigned int x_PickCritical(double prCrit) const;

    double x_TauEx(void) const {
        double tau = numeric_limits<double>::infinity();
        vector <double> mu(m_NumStates, 0);
        vector <double> sigma(m_NumStates, 0);

        for (TTransList::const_iterator j = m_TransByCat[eNormal].begin();
             j != m_TransByCat[eNormal].end();  ++j) {
            for (unsigned int i = 0;  i < m_Nu[*j].size();  ++i) {
                const SChange &c = m_Nu[*j][i];
                mu[c.m_State] += c.m_Mag * m_Rates[*j];
                sigma[c.m_State] += c.m_Mag * c.m_Mag * m_Rates[*j];
            }
        }
        //cerr << "-=| mu:";
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            //cerr << "\t" << mu[i];
            double val = max(m_Epsilon * m_X[i] / m_RateChangeBound[i],
                             1.)/fabs(mu[i]);
            //cerr << "/" << val;
            if (val < tau) {
                tau = val;
                if (tau < 0) {
                    throwError("tried to select tau < 0; most likely means "
                               "your rate function generated a negative rate");
                }
            }
            val = pow(max(m_Epsilon * m_X[i] / m_RateChangeBound[i],
                          1.),2) / sigma[i];
            if (val < tau) {
                tau = val;
                if (tau < 0) {
                    throwError("tried to select tau < 0; most likely means "
                               "your rate function generated a negative rate");
                }
            }
        }
        //cerr << endl;

        return tau;
    }

    double x_TauIm(void) const {
        if (!m_RateJacobianFunc) {
            return 0;
        }
        vector<bool> equil(m_TransCats.size(), false);
        for (TBalancedPairs::const_iterator i = m_BalancedPairs.begin();
             i != m_BalancedPairs.end();  ++i) {
            if (fabs(m_Rates[i->first] - m_Rates[i->second]) <=
                m_Delta * min(m_Rates[i->first], m_Rates[i->second])) {
                equil[i->first] = true;
                equil[i->second] = true;
            }
        }

        vector<double> mu(m_NumStates, 0);
        vector<double> sigma(m_NumStates, 0);
        for (TTransList::const_iterator j = m_TransByCat[eNormal].begin();
             j != m_TransByCat[eNormal].end();  ++j) {
            if (!equil[*j]) {
                for (unsigned int i = 0;  i < m_Nu[*j].size();  ++i) {
                    const SChange &c = m_Nu[*j][i];
                    mu[c.m_State] += c.m_Mag * m_Rates[*j];
                    sigma[c.m_State] += c.m_Mag * c.m_Mag * m_Rates[*j];
                }
            }
        }

        double tau = numeric_limits<double>::infinity();
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            double val = max(m_Epsilon * m_X[i] / m_RateChangeBound[i],
                             1.)/fabs(mu[i]);
            if (val < tau) {
                tau = val;
            }
            val = pow(max(m_Epsilon * m_X[i] / m_RateChangeBound[i],
                          1.),2) / sigma[i];
            if (val < tau) {
                tau = val;
            }
        }
        
        return tau;
    }

private:
    bool m_ExtraChecks; //turns on extra checks on rates returned by
                        //user-supplied rate function. Slower, but if
                        //the rate function does have a bug, this will
                        //give a more meaningful error message.
    int m_VerboseTracing; //trace algorithm verbosely
    bool m_RecordTransitionTimeSeries; //whether to record transitions
                                       //as a function of time in
                                       //addition to state.

    // parameters to tau leaping algorithm
    unsigned int m_Ncritical;
    double m_Nstiff;
    double m_Epsilon;
    double m_ExactThreshold;
    double m_Delta;
    unsigned int m_NumExactSteps[3];
    double m_ITLConvergenceTol;
    double m_MaxTau;
    unsigned int m_MaxSteps;

    // time-dependent variables
    double *m_T;    // *current* time
    TStates m_X;    // *current* state variables
    TRates m_Rates; // *current* rates (must be updated if m_X changes!)
    EStepType m_PrevStepType; // type of last step
    int m_LastTransition; // id of transition taken if critical/exact; -1 o.w.
    TExecutedTransitions m_ExecutedTransitions; // counts of transitions taken in current step

    // constant variables
    unsigned int m_NumStates; //total number of states
    SEXP m_VarNames;          //variable names (if any)
    TTransitions m_Nu;        //state changes caused by transitions
    TTransCats  m_TransCats;  //i.e. normal, deterministic, halting
    TTransList  m_TransByCat[4];//i.e. critical, normal, deterministic, halting
    TBalancedPairs m_BalancedPairs;
    TBools m_RealValuedVariables;
    SEXP m_RateFunc;         //R function to calculate rates as f(m_X)
    SEXP m_RateJacobianFunc; //R function to calculate Jacobian of rates as f(m_X) [optional!]
    double *m_RateChangeBound; //see Cao (2006) for details
    SEXP m_MaxTauFunc; //R function to calculate maximum leap given curr. state

    CTimeSeries m_TimeSeries;
    TTransitionTimeSeries m_TransitionTimeSeries; // record transitions in previous time interval
};


/*---------------------------------------------------------------------------*/
// PRE : m_Nu initialized
// POST: all balanced pairs of transitions identified & saved
void CStochasticEqns::x_IdentifyBalancedPairs(void) {
    for (unsigned int j1 = 0;  j1 < m_Nu.size();  ++j1) {
        for (unsigned int j2 = j1 + 1;  j2 < m_Nu.size();  ++j2) {
            if (m_Nu[j1].size() != m_Nu[j2].size()) {
                continue;
            }
            unsigned int i;
            for (i = 0;  i < m_Nu[j1].size()  &&
                     m_Nu[j1][i].m_State == m_Nu[j2][i].m_State  &&
                     m_Nu[j1][i].m_Mag == -m_Nu[j2][i].m_Mag;  ++i);
            if (i == m_Nu[j1].size()) {
                m_BalancedPairs.push_back(TBalancedPairs::value_type(j1, j2));
                if (debug) {
                    cerr << "balanced pair " << j1 << " and " << j2 << endl;
                }
            }
        }
    }
}

/*---------------------------------------------------------------------------*/
// PRE : boolean vector flagging transitions as category "cat"
// POST: appropriate transCats set
void CStochasticEqns::x_SetCat(SEXP trans, ETransCat cat) {
    if (!trans  || Rf_isNull(trans)) { return; } //NULL may be passed as a flag
    if (Rf_isLogical(trans)) {
        CRVector<bool> logic(trans);
        if (logic.size() > m_TransCats.size()) {
            throwError("length of logical vector specifying deterministic or "
                       "halting transitions is greater than the total number "
                       "of transitions!");
        }
        for (unsigned int i = 0;  i < logic.size();  ++i) {
            if (logic[i]) {
                m_TransCats[i] = cat;
                m_TransByCat[cat].push_back(i);
            }
        }
    } else {
        CRVector<int> w(PROTECT(Rf_coerceVector(trans, INTSXP)));
        UNPROTECT(1);
        for (unsigned int i = 0;  i < w.size();  ++i) {
            if (w[i] > (int) m_TransCats.size()) {
                throwError("one of your list(s) of transitions references a "
                           "transition that doesn't exist (" << w[i] << ") "
                           "when last transition is " << m_TransCats.size() <<
                           ")")
            }
            m_TransCats[w[i]-1] = cat;
            m_TransByCat[cat].push_back(w[i]-1);
        }
    }
}

/*---------------------------------------------------------------------------*/
// PRE : m_Nu initialized, deterministic transition set (if any)
// POST: all variables identified will take real values
// (i.e. either non-integer nu or modified by a deterministic transition)
void CStochasticEqns::x_IdentifyRealValuedVariables(void) {
    m_RealValuedVariables.clear();
    m_RealValuedVariables.resize(m_NumStates, false);

    for (TTransList::const_iterator j = m_TransByCat[eDeterministic].begin();
         j != m_TransByCat[eDeterministic].end();  ++j) {
        for (unsigned int i = 0;  i < m_Nu[*j].size();  ++i) {
            m_RealValuedVariables[m_Nu[*j][i].m_State] = true;
        }
    }
                    
}

/*---------------------------------------------------------------------------*/
// PRE : list of critical transitions & their total rate
// POST: one picked according to probability
unsigned int CStochasticEqns::x_PickCritical(double critRate) const {
    double r = runif(0,1);
    double d = 0;
    TTransList::const_iterator j = m_TransByCat[eCritical].begin();
    while (j != m_TransByCat[eCritical].end()) {
        d += m_Rates[*j]/critRate;
        if (d > r) {
            break;
        }
        ++j;
    }
    if (!(d >= r)) { throwError("logic error at line " << __LINE__) }
    return *j;
}

/*---------------------------------------------------------------------------*/
// PRE : time period to step; whether to clamp variables at 0
// POST: all determinisitic transitions updated by the expected amount
// (i.e. Euler method); if clamping, then negative variables set to 0.
void CStochasticEqns::x_AdvanceDeterministic(double deltaT, bool clamp) {
    for (TTransList::const_iterator j = m_TransByCat[eDeterministic].begin();
         j != m_TransByCat[eDeterministic].end();  ++j) {
        for (unsigned int i = 0;  i < m_Nu[*j].size();  ++i) {
            m_X[m_Nu[*j][i].m_State] += m_Nu[*j][i].m_Mag * m_Rates[*j] *
                deltaT;
            //clamp at zero if specified
            if (clamp  &&  m_X[m_Nu[*j][i].m_State] < 0) {
                m_X[m_Nu[*j][i].m_State] = 0;
            }
        }
        if (m_RecordTransitionTimeSeries) {
            m_ExecutedTransitions[*j] += m_Rates[*j];
        }
    }
}

/*---------------------------------------------------------------------------*/
// PRE : simulation end time; **transition rates already updated**
// POST: id of transition taken (if none, then -1) & time series updated.
void CStochasticEqns::x_SingleStepExact(double tf) {
    m_LastTransition = -1;
    m_ExecutedTransitions.assign(m_Nu.size(), 0);
    double stochRate = 0;
    double detRate = 0;
    for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
        if (m_TransCats[j] != eDeterministic) {
            stochRate += m_Rates[j];
        } else {
            detRate += m_Rates[j];
        }
    }

    double tau = stochRate > 0 ? rexp(1./stochRate) :
        detRate > 0 ? 1./detRate : tf - *m_T;
    if (stochRate == 0  ||  tau > tf - *m_T) {
        tau = tf - *m_T; // step is off end so just advance time
    } else {
        double r = runif(0,1);
        double d = 0;
        unsigned int j = 0;
        for (;  j < m_Nu.size()  &&  d < r;  ++j) {
            if (m_TransCats[j] != eDeterministic) {
                d += m_Rates[j]/stochRate;
            }
        }
        if (!(d >= r)) { throwError("logic error at line " << __LINE__) }
        --j;

        //take transition "j"
        if (m_VerboseTracing >= 1) {
            REprintf("%f: taking transition #%i\n", *m_T, j+1);
        }
        for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
            m_X[m_Nu[j][i].m_State] += m_Nu[j][i].m_Mag;
        }
        if (m_RecordTransitionTimeSeries) {
            ++m_ExecutedTransitions[j];
        }
        m_LastTransition = j;
    }

    //clamp deterministic at 0, assuming that it is unreasonable to
    //take a smaller step then exact.
    x_AdvanceDeterministic(tau, true);
    *m_T += tau;
    m_TimeSeries.push_back(STimePoint(*m_T, m_X, m_NumStates));
    if (m_RecordTransitionTimeSeries) {
        m_TransitionTimeSeries.push_back(m_ExecutedTransitions);
    }
}

/*---------------------------------------------------------------------------*/
// PRE : tau value to use for step, list of "critical" transitions
// POST: IMPLICIT tau step taken (m_X updated if so) (or overflow
// error thrown if tau was too big)
// NOTE: See equation (7) in Cao et al. (2007)
void CStochasticEqns::x_SingleStepITL(double tau) {
    if (m_VerboseTracing >= 1) {
        REprintf("%f: taking implicit step of tau = %f\n", *m_T, tau);
    }
    if (!m_RateJacobianFunc) { throwError("logic error at line " << __LINE__) }
    double *origX = new double[m_NumStates];
    double *origRates = new double[m_NumStates];
    memcpy(origX, m_X, sizeof(double)*m_NumStates);
    memcpy(origRates, m_Rates, sizeof(double)*m_NumStates);

    if (debug) {
        cerr << " origX: ";
        for (unsigned int i =0; i < m_NumStates;  ++i) {
            cerr << origX[i] << "\t";
        }
        cerr << endl;
    }

    // draw (stochastic) number of times each transition will occur
    vector<int> numTransitions(m_Nu.size(), 0);
    
    for (TTransList::const_iterator j = m_TransByCat[eNormal].begin();
         j != m_TransByCat[eNormal].end();  ++j) {
        if (m_Rates[*j]*tau > 1e8) {
            //for high rate, use normal to approx poisson.
            //should basically never yield negative, but just to
            //be sure, cap at 0
            numTransitions[*j] = max(0.,floor(rnorm(m_Rates[*j]*tau,
                                                    sqrt(m_Rates[*j]*tau))));
        } else {
            numTransitions[*j] = rpois(m_Rates[*j]*tau);
            //cerr << "nt[" << *j << "] " << numTransitions[*j] << endl;
        }
    }

    // Calculate equation (7) terms not involving x[t+tau] and call this alpha:
    //   alpha = x + nu.(P - tau/2 R(x))
    // Also initialize iterative search for x[t+tau] at expectation (reset m_X)
    double* alpha = new double[m_NumStates];
    memcpy(alpha, m_X, sizeof(double)*m_NumStates);
    for (TTransList::const_iterator j = m_TransByCat[eNormal].begin();
         j != m_TransByCat[eNormal].end();  ++j) {
        for (unsigned int k = 0;  k < m_Nu[*j].size();  ++k) {
            alpha[m_Nu[*j][k].m_State] += m_Nu[*j][k].m_Mag * 
                (numTransitions[*j] - (tau/2)*m_Rates[*j]);
            //reset m_X to expectation as our initial guess
            m_X[m_Nu[*j][k].m_State] += m_Nu[*j][k].m_Mag *
                (tau/2)*m_Rates[*j];
        }
    }
    //expectations may send states negative; clamp!
    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
        if (m_X[i] < 0) {
            m_X[i] = 0;
        }
    }

    if (debug) {
        cerr << " alpha:";
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            cerr << " " << alpha[i];
        }
        cerr << endl;
        cerr << "it " << 0 << " newX: ";
        for (unsigned int i =0; i < m_NumStates;  ++i) {
            cerr << m_X[i] << "\t";
        }
        cerr << endl;
    }

    //a few variables needed by LAPACK
    int N = m_NumStates;
    int nrhs = 1;
    int *ipiv = new int[m_NumStates];
    int info;
    double *matrixA = new double[m_NumStates*m_NumStates];
    double *matrixB = new double[m_NumStates];

    
    //Use Newton's method to solve implicit equation:
    //  Let Y = x[t+tau]
    //  F(Y) = Y - alpha - nu.((tau/2)*R(Y))
    //Solve Jacobian(F(Y0)) Y1 = -F(Y0) for Y1 to iteratively approach solution
    //This eqn expands to (I - nu.((tau/2)Jacobian(R(Y0)))) Y1 = -F(Y0) where
    //the Jacobian of rates is supplied by the user.  The term to the
    //left of Y1 is called matrix A by LAPACK and the term on right is
    //matrix B.
    //
    //Perhaps should adjust max # of iterations..
    bool converged = false;
    unsigned int c = 0;
    while (++c <= 20  &&  !converged) {
        // Check to make sure we haven't taken too big a step --
        // i.e. no state variables should go negative
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            if (m_X[i] < 0) {
                delete[] origRates;
                delete[] alpha;
                delete[] ipiv;
                delete[] matrixA;
                delete[] matrixB;
                memcpy(m_X, origX, sizeof(double)*m_NumStates);
                delete[] origX;
                throw overflow_error("tau too big");
            }
        }

        // define matrix A
        double* rateJacobian = x_CalcJacobian();
        memset(matrixA, 0, m_NumStates*m_NumStates*sizeof(double));
        for (TTransList::const_iterator j =
                 m_TransByCat[eNormal].begin();
             j != m_TransByCat[eNormal].end();  ++j) {
            for (unsigned int k = 0;  k < m_Nu[*j].size();  ++k) {
                for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                    //R stores matrices column-wise
                    //LAPACK stores matrices row-wise
                    matrixA[i*m_NumStates + m_Nu[*j][k].m_State] +=
                        m_Nu[*j][k].m_Mag * rateJacobian[(*j)*m_NumStates + i];
                }
            }
        }
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            for (unsigned int i2 = 0;  i2 < m_NumStates;  ++i2) {
                matrixA[i*m_NumStates + i2] *= -tau/2;
            }
            matrixA[i*m_NumStates + i] += 1;
        }

        // define matrix B
        // m_X is now our proposed x[t+tau].  Note that m_X has changed
        // even in our first iteration (initialized to expected value).
        x_UpdateRates();
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            matrixB[i] = alpha[i] - m_X[i];
        }
        for (TTransList::const_iterator j =
                 m_TransByCat[eNormal].begin();
             j != m_TransByCat[eNormal].end();  ++j) {
            for (unsigned int k = 0;  k < m_Nu[*j].size();  ++k) {
                matrixB[m_Nu[*j][k].m_State] += 
                    m_Nu[*j][k].m_Mag * (tau/2) * m_Rates[*j];
            }
        }


    if (debug) {
        cerr << "A:" << endl;
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            for (unsigned int i2 = 0;  i2 < m_NumStates;  ++i2) {
                cerr << matrixA[i2*m_NumStates + i] << "\t";
            }
            cerr << endl;
        }

        cerr << "B:" << endl;
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            cerr << matrixB[i] << "\t";
        }
        cerr << endl;

        cerr << "a:" << endl;
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            cerr << m_Rates[j] << "\t";
        }
        cerr << endl;
    }

        //solve linear eqn
        F77_NAME(dgesv)(&N, &nrhs, matrixA, &N, ipiv, matrixB, &N, &info); 
        if (info != 0) {
            Rf_warning("warning: lapack ran into trouble solving implicit equation");
            break;
        }
        //matrixB now contains solution (change in X)
        double normDelta = 0, normX = 0;
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            m_X[i] += matrixB[i];
            normDelta += matrixB[i]*matrixB[i];
            normX += m_X[i] * m_X[i];
        }
        //cerr << "\tNorms: " << normDelta << "\t" << normX << endl;
        converged = (normDelta < normX * m_ITLConvergenceTol);
        if (debug) {
            /*
            cerr << "Delta: ";
            for (unsigned int i =0; i < m_NumStates;  ++i) {
                cerr << matrixB[i] << "\t";
            }
            cerr << endl;
            */
            cerr << "it " << c << " newX: ";
            for (unsigned int i =0; i < m_NumStates;  ++i) {
                cerr << m_X[i] << "\t";
            }
            cerr << endl;
            /*
            x_UpdateRates();
            double t[m_NumStates];
            memcpy(t, alpha, sizeof(double)*m_NumStates);
            for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
                if (m_TransCats[j] == eNoncritical) {
                    for (unsigned int k = 0;  k < m_Nu[j].size();  ++k) {
                        t[m_Nu[j][k].m_State] +=
                            m_Nu[j][k].m_Mag * (tau/2) * m_Rates[j];
                    }
                }
            }
            cerr << "     newX: ";
            for (unsigned int i =0; i < m_NumStates;  ++i) {
                cerr << t[i] << "\t";
            }
            cerr << endl;
            */
        }

    } // end of iterating for Newton's method
    if (!converged) {
        Rf_warning("ITL solution did not converge!");
    }
    if (m_RecordTransitionTimeSeries) {
        for (TTransList::const_iterator j = m_TransByCat[eNormal].begin();
             j != m_TransByCat[eNormal].end();  ++j) {
            //round to integer # of transitions
            ++m_ExecutedTransitions[*j] +=
                floor(m_Rates[*j] + numTransitions[*j] + origRates[*j] + 0.5);
        }
    }

    //restore original rates to execute deterministic transitions
    memcpy(m_Rates, origRates, sizeof(double)*m_NumStates);
    x_AdvanceDeterministic(tau);

    delete[] origRates;
    delete[] alpha;
    delete[] ipiv;
    delete[] matrixA;
    delete[] matrixB;

    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
        if (m_X[i] < 0) {
            memcpy(m_X, origX, sizeof(double)*m_NumStates);
            delete[] origX;
            throw overflow_error("tau too big");
        }
        if (!m_RealValuedVariables[i]) {
            m_X[i] = floor(m_X[i] + 0.5); //i.e., round
        }
    }
    delete[] origX;
    *m_T += tau;
}

/*---------------------------------------------------------------------------*/
// PRE : tau value to use for step, list of "critical" transitions
// POST: EXPLICIT tau step taken (m_X updated if so) (or overflow
// error thrown if tau was too big)
void CStochasticEqns::x_SingleStepETL(double tau) {
    if (m_VerboseTracing >= 1) {
        REprintf("%f: taking explicit step of tau = %f\n", *m_T, tau);
    }
    if (m_VerboseTracing >= 2) {
        REprintf("%f:    ", *m_T);
    }
    double *origX = new double[m_NumStates];
    memcpy(origX, m_X, sizeof(double)*m_NumStates);
    for (TTransList::const_iterator j = m_TransByCat[eNormal].begin();
         j != m_TransByCat[eNormal].end();  ++j) {
        double k;
        if (m_Rates[*j]*tau > 1e8) {
            //for high rate, use normal to approx poisson.
            //should basically never yield negative, but just to
            //be sure, bound at 0
            k = max(0.,floor(rnorm(m_Rates[*j]*tau, sqrt(m_Rates[*j]*tau))));
        } else {
            k = rpois(m_Rates[*j]*tau);
        }
        if (k > 0) {
            if (m_VerboseTracing >= 2) {
                REprintf("%fx#%i ", k, *j);
            }
            for (unsigned int i = 0;  i < m_Nu[*j].size();  ++i) {
                m_X[m_Nu[*j][i].m_State] +=  k * m_Nu[*j][i].m_Mag;
            }
            if (m_RecordTransitionTimeSeries) {
                m_ExecutedTransitions[*j] += k;
            }
        }
    }
    if (m_VerboseTracing >= 2) {
        REprintf("\n");
    }
    x_AdvanceDeterministic(tau);

    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
        if (m_X[i] < 0) {
            memcpy(m_X, origX, sizeof(double)*m_NumStates);
            delete[] origX;
            throw overflow_error("tau too big");
        }
    }

    *m_T += tau;
    delete[] origX;
}

/*---------------------------------------------------------------------------*/
// PRE : time at which to end simulation; **transition rates already updated**
// POST: single adaptive tau leaping step taken & time series updated.
// Implemented from Cao Y, Gillespie DT, Petzold LR. The Journal of Chemical
// Physics (2007).
void CStochasticEqns::x_SingleStepATL(double tf) {
    m_LastTransition = -1;
    m_ExecutedTransitions.assign(m_Nu.size(), 0);
    EStepType stepType;

    //identify "critical" transitions
    double criticalRate = 0;
    double noncritRate = 0;
    {
        for (TTransList::const_iterator j =
                 m_TransByCat[eDeterministic].begin();
             j != m_TransByCat[eDeterministic].end();  ++j) {
            noncritRate += m_Rates[*j];
        }
        for (TTransList::const_iterator j =
                 m_TransByCat[eHalting].begin();
             j != m_TransByCat[eHalting].end();  ++j) {
            criticalRate += m_Rates[*j];
        }

        //reset (lop off) all non-halting criticals
        m_TransByCat[eCritical].resize(m_TransByCat[eHalting].size());
        m_TransByCat[eNormal].clear();
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] != eNormal) {
                continue;
            }
            unsigned int minTimes = numeric_limits<unsigned int>::max();
            for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                if (m_Nu[j][i].m_Mag < 0  &&
                    m_X[m_Nu[j][i].m_State]/abs(m_Nu[j][i].m_Mag) < minTimes) {
                    minTimes = m_X[m_Nu[j][i].m_State]/abs(m_Nu[j][i].m_Mag);
                }
            }
            if (minTimes < m_Ncritical) {
                criticalRate += m_Rates[j];
                m_TransByCat[eCritical].push_back(j);
            } else {
                noncritRate += m_Rates[j];
                m_TransByCat[eNormal].push_back(j);
            }
        }
    }

    if (debug) {
        cerr << "critical rate: " << criticalRate << "\t" << "noncrit rate: " << noncritRate << endl;
    }
    if (criticalRate + noncritRate == 0  &&  isinf(m_MaxTau)) {
        *m_T = tf;
        m_TimeSeries.push_back(STimePoint(*m_T, m_X, m_NumStates));
        if (m_RecordTransitionTimeSeries) {
            m_TransitionTimeSeries.push_back(m_ExecutedTransitions);
        }
        return;
    }
    if (!R_finite(criticalRate + noncritRate)) {
        throwEarlyExit("Infinite transition rate at time " << *m_T);
    }

    // calc explicit & implicit taus
    double tau1, tau2;
    double tauEx = x_TauEx();
    double tauIm = x_TauIm();
    if (debug) {
        cerr << "tauEx: " << tauEx << "  tauIm:" << tauIm << endl;
    }
    if (tauEx*m_Nstiff < tauIm) {
        stepType = eImplicit;
        tau1 = tauIm;
    } else {
        stepType = eExplicit;
        tau1 = tauEx;
    }
    if (tau1 > tf - *m_T) { //cap at the final simulation time
        tau1 = tf - *m_T;
    }
    double currMaxTau = m_MaxTau;
    if (tau1 > currMaxTau) {
        currMaxTau = m_MaxTauFunc ? min(tau1, x_CalcUserMaxTau()) : m_MaxTau;
        tau1 = currMaxTau;
        if (debug) {
            cerr << "maxtau: " << tau1 << endl;
        }
    }

    bool tauTooBig;
    do {
        tauTooBig = false;
        if (!(tau1 > 0)) { throwError("logic error at line " << __LINE__) }
        if (tau1 < m_ExactThreshold / (criticalRate + noncritRate)) {
            if (debug) {
                cerr << "Taking exact steps.. (tau1 = " << tau1 << ")" << endl;
            }
            stepType = eExact;
            for (unsigned int i = 0;
                 i < m_NumExactSteps[m_PrevStepType]  &&  *m_T < tf;  ++i) {
                if (i > 0) {
                    x_UpdateRates();
                }
                x_SingleStepExact(min(tf, *m_T + currMaxTau));
                if (m_VerboseTracing >= 2) {
                    REprintf("%f -- ", *m_T);
                    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                        REprintf("%f ", m_X[i]);
                    }
                    REprintf("\n");
                }
                if (m_LastTransition >= 0  &&
                    m_TransCats[m_LastTransition] == eHalting) {
                    return;
                }
            }
        } else {
            try { //catch exception if tauTooBig
                tau2 = (criticalRate == 0) ? numeric_limits<double>::infinity():
                    rexp(1./criticalRate);
                if (stepType == eExplicit  ||
                    (tau1 > tau2  &&  stepType == eImplicit && tau2 <= tauEx)) {
                    if (debug) {
                        cerr << "going explicit w/ tau = " << min(tau1, tau2)
                             << endl;
                    }
                    x_SingleStepETL(min(tau1, tau2));
                } else {
                    if (debug) {
                        cerr << "going implicit w/ tau = " << tau1 << endl;
                    }
                    x_SingleStepITL(tau1);
                }
                if (tau1 > tau2) { //pick one critical transition
                    unsigned int j = x_PickCritical(criticalRate);
                    m_LastTransition = j;
                    if (debug) {
                        cerr << "hittin' the critical (" << j << ")" << endl;
                    }
                    if (m_VerboseTracing >= 1) {
                        REprintf("%f:    executing critical transition #%i\n",
                                 *m_T, j+1);
                    }
                    if (m_RecordTransitionTimeSeries) {
                        ++m_ExecutedTransitions[j];
                    }
                    for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                        m_X[m_Nu[j][i].m_State] +=  m_Nu[j][i].m_Mag;
                        if (m_X[m_Nu[j][i].m_State] < 0) {
                            throwError("variable " << m_Nu[j][i].m_State+1 <<
                                       " went negative after executing "
                                       "transition " << j+1 << ".  Most likely "
                                       "either your rate calculation or "
                                       "transition matrix is flawed.");
                        }
                    }
                }

                m_TimeSeries.push_back(STimePoint(*m_T, m_X, m_NumStates));
                if (m_RecordTransitionTimeSeries) {
                    m_TransitionTimeSeries.push_back(m_ExecutedTransitions);
                }
                if (m_VerboseTracing >= 2) {
                    REprintf("%f -- ", *m_T);
                    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                        REprintf("%f ", m_X[i]);
                    }
                    REprintf("\n");
                }
            } catch (overflow_error&) { //i.e. tauTooBig exception
                tauTooBig = true;
                if (m_VerboseTracing >= 1) {
                    REprintf("%f:    tau too big; cutting in half\n", *m_T);
                }
                if (m_RecordTransitionTimeSeries) {
                    m_ExecutedTransitions.assign(m_Nu.size(), 0);
                }
                tau1 /= 2;
            }
        }
    } while (tauTooBig);

    m_PrevStepType = stepType;
}

/*---------------------------------------------------------------------------*/
// Exported C entrypoints for calling from R

extern "C" {
    SEXP simAdaptiveTau(SEXP s_x0, SEXP s_nu, SEXP s_f, SEXP s_fJacob,
                        SEXP s_params, SEXP s_tf,
                        SEXP s_deterministic, SEXP s_halting,
                        SEXP s_changebound,
                        SEXP s_tlparams, SEXP s_fMaxtau,
                        SEXP s_reportTransitions) {
        try{
        if (!Rf_isVector(s_x0)  ||  !(Rf_isReal(s_x0)  ||  Rf_isInteger(s_x0))) {
            Rf_error("invalid vector of initial values");
        }
        if (!Rf_isVectorList(s_nu)  &&
            (!Rf_isMatrix(s_nu)  ||
             INTEGER(Rf_getAttrib(s_nu, R_DimSymbol))[0] != Rf_length(s_x0))) {
            Rf_error("invalid transition specification");
        }
        if (!Rf_isFunction(s_f)) {
            Rf_error("invalid rate function");
        }
        if (!Rf_isNull(s_fJacob)  &&  !Rf_isFunction(s_fJacob)) {
            Rf_error("invalid Jacobian function");
        }
        if (!(Rf_isReal(s_tf)  ||  Rf_isInteger(s_tf))  ||  Rf_length(s_tf) != 1) {
            Rf_error("invalid final time");
        }
        if (!Rf_isVector(s_changebound)  ||  !Rf_isReal(s_changebound)  ||
            Rf_length(s_changebound) != Rf_length(s_x0)) {
            Rf_error("invalid relratechange");
        }
        if (!Rf_isNull(s_tlparams)  &&  !Rf_isVector(s_tlparams)) {
            Rf_error("tl.params must be a list");
        }
        if (!Rf_isNull(s_fMaxtau)  &&  !Rf_isFunction(s_fMaxtau)) {
            Rf_error("invalid maxTau function");
        }
        if (Rf_isNull(s_reportTransitions)  ||  !Rf_isLogical(s_reportTransitions)) {
            Rf_error("invalid value for reportTransitions");
        }

        CStochasticEqns eqns(s_x0, s_nu,
                             s_f, s_fJacob, s_params, REAL(s_changebound),
                             s_fMaxtau, s_deterministic, s_halting, s_reportTransitions);
        if (!Rf_isNull(s_tlparams)) {
            eqns.SetTLParams(s_tlparams);
        }
        try {
            eqns.EvaluateATLUntil(REAL(Rf_coerceVector(s_tf, REALSXP))[0]);
        } catch (CEarlyExit &e) {
            Rf_warning("%s", e.what());
        }
        return eqns.GetResult();
        } catch (exception &e) {
            Rf_error("%s", e.what());
            return R_NilValue;
        }
    }
    
    //-----------------------------------------------------------------------

    SEXP simExact(SEXP s_x0, SEXP s_nu, SEXP s_f, SEXP s_params, SEXP s_tf,
                  SEXP s_reportTransitions) {
        try {
        if (!Rf_isVector(s_x0)  ||  !(Rf_isReal(s_x0)  ||  Rf_isInteger(s_x0))) {
            Rf_error("invalid vector of initial values");
        }
        if (!Rf_isVectorList(s_nu)  &&
            (!Rf_isMatrix(s_nu)  ||
             INTEGER(Rf_getAttrib(s_nu, R_DimSymbol))[0] != Rf_length(s_x0))) {
            Rf_error("invalid transition specification");
        }
        if (!Rf_isFunction(s_f)) {
            Rf_error("invalid rate function");
        }
        if (!(Rf_isReal(s_tf)  ||  Rf_isInteger(s_tf))  ||  Rf_length(s_tf) != 1) {
            Rf_error("invalid final time");
        }
        if (Rf_isNull(s_reportTransitions)  ||  !Rf_isLogical(s_reportTransitions)) {
            Rf_error("invalid value for reportTransitions");
        }

        CStochasticEqns eqns(s_x0, s_nu,
                             s_f, NULL, s_params, NULL, NULL,
                             R_NilValue, R_NilValue, s_reportTransitions);
        try {
            eqns.EvaluateExactUntil(REAL(Rf_coerceVector(s_tf, REALSXP))[0]);
        } catch (CEarlyExit &e) {
            Rf_warning("%s", e.what());
        }
        return eqns.GetResult();
        } catch (exception &e) {
            Rf_error("%s", e.what());
            return R_NilValue;
        }
    }

    const R_CallMethodDef callMethods[] = {
	{"simAdaptiveTau", (DL_FUNC)&simAdaptiveTau, 12},
	{"simExact", (DL_FUNC)&simExact, 6},
	{NULL, NULL, 0}
    };
    void R_init_adaptivetau(DllInfo *dll) {
	R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
    }

}
