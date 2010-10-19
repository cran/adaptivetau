/* $Id: adaptivetau.cpp 162 2010-11-12 17:28:37Z pjohnson $
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
    --------------------------------------------------------------------------
*/

#include <iostream>
#include <vector>
#include <limits>
#include <cassert>
#include <cmath>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>

using namespace std;

enum EStepType {
    eExact = 0,
    eExplicit,
    eImplicit
};

bool debug = false;

class CStochasticEqns {
public:
    CStochasticEqns(SEXP initVal, int *nu, unsigned int numTrans,
                    SEXP rateFunc, SEXP rateJacobianFunc,
                    SEXP params, double* changeBound, SEXP maxTauFunc) {
        // copy initial values into new vector (keeping in SEXP vector
        // allows easy calling of R function to calculate rates)
        m_NumStates = GET_LENGTH(initVal);
        SEXP x;
        PROTECT(x = allocVector(REALSXP, m_NumStates));
        copyVector(x, initVal);
        if (isNull(GET_NAMES(initVal))) {
            m_VarNames = NULL;
        } else {
            SEXP namesO = GET_NAMES(initVal);
            PROTECT(m_VarNames = allocVector(VECSXP, GET_LENGTH(namesO)));
            copyVector(m_VarNames, namesO);
            SET_NAMES(x, m_VarNames);
        }
        m_X = REAL(x);

        // copy full-size Nu matrix into sparse matrix
        m_Nu.resize(numTrans);
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            for (unsigned int j = 0;  j < numTrans;  ++j) {
                if (nu[j*m_NumStates + i] != 0) {
                    m_Nu[j].push_back(SChange(i, nu[j*m_NumStates + i]));
                }
            }
        }
        m_TransCats.resize(numTrans, eNoncritical);

        // needed for ITL
        x_IdentifyBalancedPairs();
        x_IdentifyRealValuedVariables();

        // prepare R function for evaluation by setting up arguments
        // (current X values, parameters, current time)
        SEXP s_time;
        PROTECT(s_time = allocVector(REALSXP, 1));
        m_T = REAL(s_time);
        PROTECT(m_RateFunc = lang4(rateFunc, x, params, s_time));
        if (isNull(rateJacobianFunc)) {
            m_RateJacobianFunc = NULL;
        } else {
            PROTECT(m_RateJacobianFunc = lang4(rateJacobianFunc, x,
                                               params, s_time));
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
        m_ITLConvergenceTol = 0.1;
        m_MaxTau = numeric_limits<double>::infinity();
        m_MaxSteps = 0; // special case 0 == no limit

        //useful additional parameters
        m_RateChangeBound = changeBound;
        if (isNull(maxTauFunc)) {
            m_MaxTauFunc = NULL;
        } else {
            PROTECT(m_MaxTauFunc = lang4(maxTauFunc, x, params,s_time));
        }
        
        *m_T = 0;
        m_PrevStepType = eExact;
        GetRNGstate();
    }
    ~CStochasticEqns(void) {
        int cnt = 3;
        if (m_RateJacobianFunc != NULL) {
            ++cnt;
        }
        if (m_Rates != NULL) {
            ++cnt;
        }
        if (m_MaxTauFunc != NULL) {
            ++cnt;
        }
        if (m_VarNames != NULL) {
            ++cnt;
        }
        UNPROTECT(cnt);
        PutRNGstate();
    }
    void SetTLParams(SEXP list) {
        SEXP names = GET_NAMES(list);

        for (int i = 0;  i < GET_LENGTH(names);  ++i) {
            if (!isReal(VECTOR_ELT(list, i))) {
                warning("ignoring non-numeric value for parameter '%s'",
                        CHAR(VECTOR_ELT(names, i)));
            } else if (strcmp("epsilon", CHAR(CHARACTER_DATA(names)[i])) == 0) {
                m_Epsilon = REAL(VECTOR_ELT(list, i))[0];
            } else if (strcmp("delta", CHAR(CHARACTER_DATA(names)[i])) == 0) {
                m_Delta = REAL(VECTOR_ELT(list, i))[0];
            } else if (strcmp("maxtau", CHAR(CHARACTER_DATA(names)[i])) == 0) {
                m_MaxTau = REAL(VECTOR_ELT(list, i))[0];
            } else {
                warning("ignoring unknown parameter '%s'",
                        CHAR(CHARACTER_DATA(names)[i]));
            }
        }
    }

    void SetDeterministic(int *det, unsigned int n) {
        if (n != m_Nu.size()  &&  n > 1) {
            warning("mismatch between length of logical vector specifying deterministic transitions and total number of transitions; ignoring deterministic!");
            return;
        }
        for (unsigned int i = 0;  i < n;  ++i) {
            if (det[i]) {
                m_TransCats[i] = eDeterministic;
            }
        }
        x_IdentifyRealValuedVariables();
    }

    void EvaluateUntil(double tF) {
        unsigned int c = 0;
        while (*m_T < tF  &&  (m_MaxSteps == 0 || ++c < m_MaxSteps)) {
            x_SingleStepATL(tF);
            R_CheckUserInterrupt();
        }
    }
    SEXP GetTimeSeriesSEXP(void) const {
        SEXP res;
        PROTECT(res = allocMatrix(REALSXP, m_TimeSeries.size(), m_NumStates+1));
        double *rvals = REAL(res);
        for (unsigned int t = 0;  t < m_TimeSeries.size();  ++t) {
            rvals[t] = m_TimeSeries[t].m_T;
            for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                rvals[(i+1) * m_TimeSeries.size() + t] = m_TimeSeries[t].m_X[i];
            }
        }

        SEXP dimnames, colnames;
        PROTECT(dimnames = allocVector(VECSXP, 2));
        PROTECT(colnames = allocVector(VECSXP, m_NumStates+1));
        SET_VECTOR_ELT(dimnames, 1, colnames);
        SET_VECTOR_ELT(colnames, 0, mkChar("time"));
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            if (m_VarNames  &&  (unsigned int)GET_LENGTH(m_VarNames) > i) {
                SET_VECTOR_ELT(colnames, i+1,
                               CHARACTER_DATA(m_VarNames)[i]);
            } else {
                char name[10];
                snprintf(name, 10, "x%i", i+1);
                SET_VECTOR_ELT(colnames, i+1, mkChar(name));
            }
        }
        SET_DIMNAMES(res, dimnames);

        UNPROTECT(3);
        return res;
    }

protected:
    enum ETransCat {
        eCritical,
        eNoncritical,
        eDeterministic
    };
    typedef vector<ETransCat> TTransCats;
    typedef vector<pair<unsigned int, unsigned int> > TBalancedPairs;
    typedef vector<bool> TBools;
    typedef double* TStates;
    typedef double* TRates;
    struct SChange {
        SChange(short int s, short int m) : m_State(s), m_Mag(m) {}
        short int m_State;
        short int m_Mag;
    };
    typedef vector< vector<SChange> > TTransitions;
    struct STimePoint {
        STimePoint(double t, double *x, int n) {
            m_T = t;
            m_X = new double[n];
            memcpy(m_X, x, n*sizeof(*x));
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

protected:
    void x_IdentifyBalancedPairs(void);
    void x_IdentifyRealValuedVariables(void);

    void x_AdvanceDeterministic(double deltaT);
    void x_SingleStepExact(double tf);
    bool x_SingleStepETL(double tau);
    bool x_SingleStepITL(double tau);
    void x_SingleStepATL(double tf);

    void x_UpdateRates(void) {
        //not sure if this protect/unprotect block is necessary, but
        //seems better to err on the safe side
        if (m_Rates != NULL) { 
            UNPROTECT(1);
        }
        SEXP res = eval(m_RateFunc, NULL);
        PROTECT(res);
        m_Rates = REAL(res);
        if ((unsigned int) GET_LENGTH(res) != m_Nu.size()) {
            error("invalid rate function -- returned number of rates is "
                  "not the same as specified by the transition matrix! "
                  "(%i versus %i)", GET_LENGTH(res), m_Nu.size());
        }
    }
    double* x_CalcJacobian(void) {
        SEXP res = eval(m_RateJacobianFunc, NULL);
        if ((unsigned int)GET_LENGTH(res) != m_Nu.size()*m_NumStates) {
            error("invalid Jacobian function -- did not return enough rates!");
        }
        return REAL(res);
    }
    double x_CalcUserMaxTau(void) {
        assert(m_MaxTauFunc);
        SEXP res = eval(m_MaxTauFunc, NULL);
        if (GET_LENGTH(res) != 1  || !isReal(res)) {
            error("invalid return value from maxTau function");
        }
        return REAL(res)[0];
    }

    unsigned int x_PickCritical(double prCrit) const;

    double x_TauEx(void) const {
        double tau = numeric_limits<double>::infinity();
        vector <double> mu(m_NumStates, 0);
        vector <double> sigma(m_NumStates, 0);

        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] != eCritical) {
                for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                    const SChange &c = m_Nu[j][i];
                    mu[c.m_State] += c.m_Mag * m_Rates[j];
                    sigma[c.m_State] += c.m_Mag * c.m_Mag * m_Rates[j];
                }
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
                    error("tried to select tau < 0; most likely means your "
                          "rate function generated a negative rate");
                }
            }
            val = pow(max(m_Epsilon * m_X[i] / m_RateChangeBound[i],
                          1.),2) / sigma[i];
            if (val < tau) {
                tau = val;
                if (tau < 0) {
                    error("tried to select tau < 0; most likely means your "
                          "rate function generated a negative rate");
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
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] != eCritical  &&  !equil[j]) {
                for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                    const SChange &c = m_Nu[j][i];
                    mu[c.m_State] += c.m_Mag * m_Rates[j];
                    sigma[c.m_State] += c.m_Mag * c.m_Mag * m_Rates[j];
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
    unsigned int m_Ncritical;
    double m_Nstiff;
    double m_Epsilon;
    double m_ExactThreshold;
    double m_Delta;
    unsigned int m_NumExactSteps[3];
    double m_ITLConvergenceTol;
    double m_MaxTau;
    unsigned int m_MaxSteps;

    TRates m_Rates; // *current* rates (must be updated if m_X changes!)
    double *m_T;    // *current* time
    TBalancedPairs m_BalancedPairs;
    TBools m_RealValuedVariables;
    EStepType m_PrevStepType;

    TStates m_X;              //current state
    unsigned int m_NumStates; //total number of states
    SEXP m_VarNames;          //variable names (if any)
    TTransitions m_Nu;        //state changes caused by transition
    TTransCats m_TransCats; //inc. whether transition is deterministic
    SEXP m_RateFunc; //R function to calculate rates as f(m_X)
    SEXP m_RateJacobianFunc; //R function to calculate Jacobian of rates as f(m_X) [optional!]
    double *m_RateChangeBound; //see Cao (2006) for details
    SEXP m_MaxTauFunc; //R function to calculate maximum leap given curr. state

    CTimeSeries m_TimeSeries;
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
// PRE : m_Nu initialized, deterministic transition set (if any)
// POST: all variables identified will take real values
// (i.e. either non-integer nu or modified by a deterministic transition)
void CStochasticEqns::x_IdentifyRealValuedVariables(void) {
    m_RealValuedVariables.clear();
    m_RealValuedVariables.resize(m_NumStates, false);

    for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
        if (m_TransCats[j] == eDeterministic) {
            for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                m_RealValuedVariables[m_Nu[j][i].m_State] = true;
            }
        } else {
            for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                if (m_Nu[j][i].m_Mag - trunc(m_Nu[j][i].m_Mag) > 1e-5) {
                    m_RealValuedVariables[m_Nu[j][i].m_State] = true;
                }
            }
        }
    }
                    
}

/*---------------------------------------------------------------------------*/
// PRE : list of critical transitions & their total rate
// POST: one picked according to probability
unsigned int CStochasticEqns::x_PickCritical(double critRate) const {
    if (critRate == 0) {
        error("someting screwy -- complain to author");
    }
    double r = runif(0,1);
    double d = 0;
    unsigned int j;
    for (j = 0;  j < m_Nu.size()  &&  d < r;  ++j) {
        if (m_TransCats[j] == eCritical) {
            d += m_Rates[j]/critRate;
        }
    }
    assert(d >= r); //otherwise logic error..
    return j-1;
}

/*---------------------------------------------------------------------------*/
// PRE : time period to step
// POST: all determinisitic transitions updated by the expected amount
void CStochasticEqns::x_AdvanceDeterministic(double deltaT) {
    for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
        if (m_TransCats[j] == eDeterministic) {
            for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                m_X[m_Nu[j][i].m_State] += m_Nu[j][i].m_Mag * m_Rates[j] *
                    deltaT;
            }
        }
    }
}

/*---------------------------------------------------------------------------*/
// PRE : simulation end time
// POST: a *single* transition taken (no approximation necessary)
void CStochasticEqns::x_SingleStepExact(double tf) {
    double stochRate = 0;
    for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
        if (m_TransCats[j] != eDeterministic) {
            stochRate += m_Rates[j];
        }
    }
    if (stochRate == 0) {
        *m_T = numeric_limits<double>::infinity();
        return;
    }
    double tau = rexp(1./stochRate);
    if (tau > tf - *m_T) {
        tau = tf - *m_T;
    } else { // only take step if don't go off end
        double r = runif(0,1);
        double d = 0;
        unsigned int j;
        for (j = 0;  j < m_Nu.size()  &&  d < r;  ++j) {
            if (m_TransCats[j] != eDeterministic) {
                d += m_Rates[j]/stochRate;
            }
        }
        assert(d >= r); //otherwise logic error..
        --j;

        //take transition "j"
        for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
            m_X[m_Nu[j][i].m_State] += m_Nu[j][i].m_Mag;
        }
    }

    x_AdvanceDeterministic(tau);
    //clamp deterministic at 0, assuming that it is unreasonable to
    //take a smaller step then exact.
    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
        if (m_X[i] < 0) {
            m_X[i] = 0;
        }
    }
    *m_T += tau;
}

/*---------------------------------------------------------------------------*/
// PRE : tau value to use for step, list of "critical" transitions
// POST: whether single IMPLICIT tau step was successfully taken (m_X
// updated if so)
bool CStochasticEqns::x_SingleStepITL(double tau) {
    assert(m_RateJacobianFunc != NULL);
    double *origX = new double[m_NumStates];
    double *origRates = new double[m_NumStates];
    memcpy(origX, m_X, sizeof(double)*m_NumStates);
    memcpy(origRates, m_Rates, sizeof(double)*m_NumStates);

    vector<int> numTransitions(m_Nu.size(), 0);
    for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
        if (m_TransCats[j] == eNoncritical) {
            if (m_Rates[j]*tau > 1e8) {
                //for high rate, use normal to approx poisson.
                //should basically never yield negative, but just to
                //be sure, cap at 0
                numTransitions[j] = max(0.,floor(rnorm(m_Rates[j]*tau, sqrt(m_Rates[j]*tau))));
            } else {
                numTransitions[j] = rpois(m_Rates[j]*tau);
            }
        }
    }

    //all terms not involving x[t+tau]
    double* alpha = new double[m_NumStates];
    memcpy(alpha, m_X, sizeof(double)*m_NumStates);
    for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
        if (m_TransCats[j] == eNoncritical) {
            for (unsigned int k = 0;  k < m_Nu[j].size();  ++k) {
                alpha[m_Nu[j][k].m_State] += m_Nu[j][k].m_Mag * 
                    (numTransitions[j] - (tau/2)*m_Rates[j]);
                //reset m_X to expectation as our initial guess
                m_X[m_Nu[j][k].m_State] += m_Nu[j][k].m_Mag *
                    (tau/2)*m_Rates[j];
            }
        }
    }
    /*
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
    */
    //a few variables needed by LAPACK
    int N = m_NumStates;
    int nrhs = 1;
    int *ipiv = new int[m_NumStates];
    int info;
    double *matrixA = new double[m_NumStates*m_NumStates];
    double *matrixB = new double[m_NumStates];

    
    //Use Newton's method to solve implicit equation
    //perhaps should adjust max # of iterations..
    bool converged = false;
    unsigned int c = 0;
    while (++c <= 20  &&  !converged) {
        //cerr << "implicit iteration #" << c << endl;
        memset(matrixA, 0, m_NumStates*m_NumStates*sizeof(double));
        double *rateJacobian;
        rateJacobian = x_CalcJacobian();
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] == eNoncritical) {
                for (unsigned int k = 0;  k < m_Nu[j].size();  ++k) {
                    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                        //must be stored as col, row for LAPACK
                        matrixA[i*m_NumStates + m_Nu[j][k].m_State] +=
                            m_Nu[j][k].m_Mag * rateJacobian[j*m_NumStates + i];
                    }
                }
            }
        }
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            for (unsigned int i2 = 0;  i2 < m_NumStates;  ++i2) {
                matrixA[i*m_NumStates + i2] *= -tau/2;
            }
            matrixA[i*m_NumStates + i] += 1;
        }

        x_UpdateRates();
        for (unsigned int i = 0;  i < m_NumStates;  ++i) {
            matrixB[i] = alpha[i] - m_X[i];
        }
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] == eNoncritical) {
                for (unsigned int k = 0;  k < m_Nu[j].size();  ++k) {
                    matrixB[m_Nu[j][k].m_State] += 
                        m_Nu[j][k].m_Mag * (tau/2) * m_Rates[j];
                }
            }
        }

        /*
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
        */
        //solve eqn
        F77_NAME(dgesv)(&N, &nrhs, matrixA, &N, ipiv, matrixB, &N, &info); 
        if (info != 0) {
            cerr << "warning: lapack ran into trouble solving implicit equation"
                 << endl;
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

    }
    if (!converged) {
        cerr << "ITL solution did not converge!" << endl;
    } else {
        x_UpdateRates();
        memcpy(m_X, alpha, sizeof(double)*m_NumStates);
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] == eNoncritical) {
                for (unsigned int k = 0;  k < m_Nu[j].size();  ++k) {
                    m_X[m_Nu[j][k].m_State] +=
                        m_Nu[j][k].m_Mag * (tau/2) * m_Rates[j];
                }
            }
        }
    }
    memcpy(m_Rates, origRates, sizeof(double)*m_NumStates);
    x_AdvanceDeterministic(tau);

    delete[] origRates;
    delete[] alpha;
    delete[] ipiv;
    delete[] matrixA;
    delete[] matrixB;

    bool tauTooBig = false;
    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
        if (m_X[i] < 0) {
            tauTooBig = true;
            break;
        }
        if (!m_RealValuedVariables[i]) {
            m_X[i] = round(m_X[i]);
        }
    }
    if (tauTooBig) {
        memcpy(m_X, origX, sizeof(double)*m_NumStates);
        delete[] origX;
        return false;
    }
    delete[] origX;
    *m_T += tau;
    return true;
}

/*---------------------------------------------------------------------------*/
// PRE : tau value to use for step, list of "critical" transitions
// POST: whether single EXPLICIT tau step was successfully taken (m_X
// updated if so)
bool CStochasticEqns::x_SingleStepETL(double tau) {
    double *origX = new double[m_NumStates];
    memcpy(origX, m_X, sizeof(double)*m_NumStates);
    for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
        if (m_TransCats[j] == eNoncritical) {
            double k;
            if (m_Rates[j]*tau > 1e8) {
                //for high rate, use normal to approx poisson.
                //should basically never yield negative, but just to
                //be sure, cap at 0
                k = max(0.,floor(rnorm(m_Rates[j]*tau, sqrt(m_Rates[j]*tau))));
            } else {
                k = rpois(m_Rates[j]*tau);
            }
            for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                m_X[m_Nu[j][i].m_State] +=  k * m_Nu[j][i].m_Mag;
            }
        } else if (m_TransCats[j] == eDeterministic) {
            for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                m_X[m_Nu[j][i].m_State] += m_Nu[j][i].m_Mag * m_Rates[j] *
                    tau;
            }
        }
    }

    bool tauTooBig = false;
    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
        if (m_X[i] < 0) {
            tauTooBig = true;
            break;
        }
    }
    if (tauTooBig) {
        memcpy(m_X, origX, sizeof(double)*m_NumStates);
        delete[] origX;
        return false;
    }

    *m_T += tau;
    delete[] origX;
    return true;
}

/*---------------------------------------------------------------------------*/
// PRE : time at which to end simulation
// POST: single adaptive tau leaping step taken & time series updated.
// Implemented from Cao Y, Gillespie DT, Petzold LR. The Journal of Chemical
// Physics (2007).
void CStochasticEqns::x_SingleStepATL(double tf) {
    //debug check
    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
        if (m_X[i] < 0) {
            error("negative variable: %i is %f", i+1, m_X[i]);
        } else if (isnan(m_X[i])) {
            error("NaN variable: %i is %f", i+1, m_X[i]);
        }
    }

    x_UpdateRates();
    EStepType stepType;

    //identify "critical" transitions
    double criticalRate = 0;
    double noncritRate = 0;
    {
        for (unsigned int j = 0;  j < m_Nu.size();  ++j) {
            if (m_TransCats[j] == eDeterministic) {
                noncritRate += m_Rates[j];
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
                m_TransCats[j] = eCritical;
                criticalRate += m_Rates[j];
            } else {
                m_TransCats[j] = eNoncritical;
                noncritRate += m_Rates[j];
            }
            if (m_Rates[j] < 0) {
                error("at time %f, rate of transition %i is %f (less than 0)!", *m_T, j+1, m_Rates[j]);
            } else if (isnan(m_Rates[j])) {
                error("at time %f, rate of transition %i is not a number (NA?)", *m_T, j+1, m_Rates[j]);
            }
        }
    }

    if (debug) {
        cerr << "critical rate: " << criticalRate << "\t" << "noncrit rate: " << noncritRate << endl;
    }
    if (criticalRate + noncritRate == 0) {
        *m_T = tf;//numeric_limits<double>::infinity();
        m_TimeSeries.push_back(STimePoint(*m_T, m_X, m_NumStates));
        return;
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
    if (tau1 > m_MaxTau) {
        tau1 = m_MaxTauFunc ? min(tau1, x_CalcUserMaxTau()) : m_MaxTau;
        if (debug) {
            cerr << "maxtau: " << tau1 << " (" <<
                (m_MaxTauFunc ? x_CalcUserMaxTau() : m_MaxTau) << ")" << endl;
        }
    }

    bool tauTooBig;
    do {
        tauTooBig = false;
        if (tau1 < 0) {
            error("uhoh.. tau1 < 0!\n");
        }
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
                x_SingleStepExact(tf);
                if (isinf(*m_T)) { //signal that rates = 0
                    *m_T = tf;
                }
                if (debug) {
                    cerr << *m_T << " -- ";
                    for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                        cerr << m_X[i] << " ";
                    }
                    cerr << endl;
                }
                m_TimeSeries.push_back(STimePoint(*m_T, m_X, m_NumStates));
            }
        } else {
            tau2 = (criticalRate == 0) ? numeric_limits<double>::infinity() :
                rexp(1./criticalRate);
            if (stepType == eExplicit  ||
                (tau1 > tau2  &&  stepType == eImplicit  &&  tau2 <= tauEx)) {
                if (debug) {
                    cerr << "going explicit w/ tau = " << min(tau1, tau2) << endl;
                }
                tauTooBig = !x_SingleStepETL(min(tau1, tau2));
            } else {
                if (1) {
                    cerr << "going implicit w/ tau = " << tau1 << endl;
                }
                tauTooBig = !x_SingleStepITL(tau1);
            }
            if (!tauTooBig) {
                if (tau1 > tau2) { //pick one critical transition
                    unsigned int j = x_PickCritical(criticalRate);
                    if (debug) {
                        cerr << "hittin' the critical (" << j << ")" << endl;
                    }
                    for (unsigned int i = 0;  i < m_Nu[j].size();  ++i) {
                        m_X[m_Nu[j][i].m_State] +=  m_Nu[j][i].m_Mag;
                        if (m_X[m_Nu[j][i].m_State] < 0) {
                            error("variable %i went negative after executing "
                                  "transition %i.  Most likely either your "
                                  "rate calculation or transition matrix is "
                                  "flawed.", m_Nu[j][i].m_State+1, j+1);
                        }
                    }
                }

                m_TimeSeries.push_back(STimePoint(*m_T, m_X, m_NumStates));
            }
            if (debug) {
                cerr << *m_T << " -- ";
                for (unsigned int i = 0;  i < m_NumStates;  ++i) {
                    cerr << m_X[i] << " ";
                }
                cerr << endl;
            }
        }
        if (tauTooBig) {
            if (debug) {
                cerr << "whoa.. knock that tau back down" << endl;
            }
            tau1 /= 2;
        }
    } while (tauTooBig);

    m_PrevStepType = stepType;
}

/*---------------------------------------------------------------------------*/
// Exported C entrypoints for calling from R

extern "C" {
    SEXP adaptiveTau(SEXP s_x0, SEXP s_nu, SEXP s_f, SEXP s_fJacob,
                     SEXP s_params, SEXP s_tf,
                     SEXP s_deterministic, SEXP s_changebound,
                     SEXP s_tlparams, SEXP s_fMaxtau) {
        if (!isVector(s_x0)  ||  !isReal(s_x0)) {
            error("invalid vector of initial values");
        }
        if (!isMatrix(s_nu)  ||  !isInteger(s_nu)) {
            error("invalid transitions matrix");
        }
        if (!isFunction(s_f)) {
            error("invalid rate function");
        }
        if (!isNull(s_fJacob)  &&  !isFunction(s_fJacob)) {
            error("invalid Jacobian function");
        }
        if (!isReal(s_tf)) {
            error("invalid final time");
        }
        if (GET_LENGTH(s_x0) != INTEGER(GET_DIM(s_nu))[0]) {
            error("mismatch between number of state variables and number of rows in transition matrix");
        }
        if (!isVector(s_deterministic)  ||  !isLogical(s_deterministic)) {
            error("invalid deterministic parameter -- must be logical vector");
        }
        if (!isVector(s_changebound)  ||  !isReal(s_changebound)  ||
            GET_LENGTH(s_changebound) != GET_LENGTH(s_x0)) {
            error("invalid relratechange");
        }
        if (!isVector(s_tlparams)) {
            error("tl.params must be a vector");
        }
        if (!isNull(s_fMaxtau)  &&  !isFunction(s_fMaxtau)) {
            error("invalid maxTau function");
        }

        CStochasticEqns eqns(s_x0, INTEGER(s_nu), INTEGER(GET_DIM(s_nu))[1],
                             s_f, s_fJacob, s_params, REAL(s_changebound),
                             s_fMaxtau);
        eqns.SetTLParams(s_tlparams);
        eqns.SetDeterministic(LOGICAL(s_deterministic),
                              GET_LENGTH(s_deterministic));
        eqns.EvaluateUntil(REAL(s_tf)[0]);
        return eqns.GetTimeSeriesSEXP();
    }
    
    R_CallMethodDef callMethods[] = {
        {"ssa.adaptivetau", (DL_FUNC)&adaptiveTau, 10},
        {NULL,NULL, 0}
    };

    void R_init_adaptivetau(DllInfo *dll)
    {
        R_registerRoutines(dll,NULL,callMethods,NULL,NULL);
    }

}
