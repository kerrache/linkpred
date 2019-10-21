/*
 * cgdescent.hpp
 *
 *  Created on: Aug 8, 2016
 *      Author: Said Kerrache
 */

#ifndef CGDESCENT_HPP_
#define CGDESCENT_HPP_

#include <limits.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include "linkpred/numerical/cg/cgblas.hpp"

namespace CG {
typedef long int INT;
const INT INT_INF = LONG_MAX;
constexpr double INF = DBL_MAX;

constexpr double ZERO = 0;
constexpr double ONE = 1;
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/*============================================================================
 cg_parameter is a structure containing parameters used in cg_descent
 cg_default assigns default values to these parameters */
typedef struct cg_parameter_struct /* user controlled parameters */
{
	/*============================================================================
	 parameters that the user may wish to modify
	 ----------------------------------------------------------------------------*/
	/* T => print final statistics
	 F => no printout of statistics */
	int PrintFinal;

	/* Level 0  = no printing), ... , Level 3 = maximum printing */
	int PrintLevel;

	/* T => print parameters values
	 F => do not display parmeter values */
	int PrintParms;

	/* T => use LBFGS
	 F => only use L-BFGS when memory >= n */
	int LBFGS;

	/* number of vectors stored in memory */
	int memory;

	/* SubCheck and SubSkip control the frequency with which the subspace
	 condition is checked. It is checked for SubCheck*mem iterations and
	 if not satisfied, then it is skipped for Subskip*mem iterations
	 and Subskip is doubled. Whenever the subspace condition is statisfied,
	 SubSkip is returned to its original value. */
	int SubCheck;
	int SubSkip;

	/* when relative distance from current gradient to subspace <= eta0,
	 enter subspace if subspace dimension = mem */
	double eta0;

	/* when relative distance from current gradient to subspace >= eta1,
	 leave subspace */
	double eta1;

	/* when relative distance from current direction to subspace <= eta2,
	 always enter subspace (invariant space) */
	double eta2;

	/* T => use approximate Wolfe line search
	 F => use ordinary Wolfe line search, switch to approximate Wolfe when
	 |f_k+1-f_k| < AWolfeFac*C_k, C_k = average size of cost  */
	int AWolfe;
	double AWolfeFac;

	/* factor in [0, 1] used to compute average cost magnitude C_k as follows:
	 Q_k = 1 + (Qdecay)Q_k-1, Q_0 = 0,  C_k = C_k-1 + (|f_k| - C_k-1)/Q_k */
	double Qdecay;

	/* terminate after nslow iterations without strict improvement in
	 either function value or gradient */
	int nslow;

	/* Stop Rules:
	 T => ||proj_grad||_infty <= max(grad_tol,initial ||grad||_infty*StopFact)
	 F => ||proj_grad||_infty <= grad_tol*(1 + |f_k|) */
	int StopRule;
	double StopFac;

	/* T => estimated error in function value is eps*Ck,
	 F => estimated error in function value is eps */
	int PertRule;
	double eps;

	/* factor by which eps grows when line search fails during contraction */
	double egrow;

	/* T => attempt quadratic interpolation in line search when
	 |f_k+1 - f_k|/f_k <= QuadCutoff
	 F => no quadratic interpolation step */
	int QuadStep;
	double QuadCutOff;

	/* maximum factor by which a quad step can reduce the step size */
	double QuadSafe;

	/* T => when possible, use a cubic step in the line search */
	int UseCubic;

	/* use cubic step when |f_k+1 - f_k|/|f_k| > CubicCutOff */
	double CubicCutOff;

	/* |f| < SmallCost*starting cost => skip QuadStep and set PertRule = FALSE*/
	double SmallCost;

	/* T => check that f_k+1 - f_k <= debugtol*C_k
	 F => no checking of function values */
	int debug;
	double debugtol;

	/* if step is nonzero, it is the initial step of the initial line search */
	double step;

	/* abort cg after maxit iterations */
	INT maxit;

	/* maximum number of times the bracketing interval grows during expansion */
	int ntries;

	/* maximum factor secant step increases stepsize in expansion phase */
	double ExpandSafe;

	/* factor by which secant step is amplified during expansion phase
	 where minimizer is bracketed */
	double SecantAmp;

	/* factor by which rho grows during expansion phase where minimizer is
	 bracketed */
	double RhoGrow;

	/* maximum number of times that eps is updated */
	int neps;

	/* maximum number of times the bracketing interval shrinks */
	int nshrink;

	/* maximum number of iterations in line search */
	int nline;

	/* conjugate gradient method restarts after (n*restart_fac) iterations */
	double restart_fac;

	/* stop when -alpha*dphi0 (estimated change in function value) <= feps*|f|*/
	double feps;

	/* after encountering nan, growth factor when searching for
	 a bracketing interval */
	double nan_rho;

	/* after encountering nan, decay factor for stepsize */
	double nan_decay;

	/*============================================================================
	 technical parameters which the user probably should not touch
	 ----------------------------------------------------------------------------*/
	double delta; /* Wolfe line search parameter */
	double sigma; /* Wolfe line search parameter */
	double gamma; /* decay factor for bracket interval width */
	double rho; /* growth factor when searching for initial
	 bracketing interval */
	double psi0; /* factor used in starting guess for iteration 1 */
	double psi_lo; /* in performing a QuadStep, we evaluate at point
	 between [psi_lo, psi_hi]*psi2*previous step */
	double psi_hi;
	double psi1; /* for approximate quadratic, use gradient at
	 psi1*psi2*previous step for initial stepsize */
	double psi2; /* when starting a new cg iteration, our initial
	 guess for the line search stepsize is
	 psi2*previous step */
	int AdaptiveBeta; /* T => choose beta adaptively, F => use theta */
	double BetaLower; /* lower bound factor for beta */
	double theta; /* parameter describing the cg_descent family */
	double qeps; /* parameter in cost error for quadratic restart
	 criterion */
	double qrule; /* parameter used to decide if cost is quadratic */
	int qrestart; /* number of iterations the function should be
	 nearly quadratic before a restart */
} cg_parameter;

typedef struct cg_stats_struct /* statistics returned to user */
{
	double f; /*function value at solution */
	double gnorm; /* max abs component of gradient */
	INT iter; /* number of iterations */
	INT IterSub; /* number of subspace iterations */
	INT NumSub; /* total number subspaces */
	INT nfunc; /* number of function evaluations */
	INT ngrad; /* number of gradient evaluations */
} cg_stats;

typedef struct cg_com_struct /* common variables */
{
	/* parameters computed by the code */
	INT n; /* problem dimension, saved for reference */
	INT nf; /* number of function evaluations */
	INT ng; /* number of gradient evaluations */
	int QuadOK; /* T (quadratic step successful) */
	int UseCubic; /* T (use cubic step) F (use secant step) */
	int neps; /* number of time eps updated */
	int PertRule; /* T => estimated error in function value is eps*Ck,
	 F => estimated error in function value is eps */
	int QuadF; /* T => function appears to be quadratic */
	double SmallCost; /* |f| <= SmallCost => set PertRule = F */
	double alpha; /* stepsize along search direction */
	double f; /* function value for step alpha */
	double df; /* function derivative for step alpha */
	double fpert; /* perturbation is eps*|f| if PertRule is T */
	double eps; /* current value of eps */
	double tol; /* computing tolerance */
	double f0; /* old function value */
	double df0; /* old derivative */
	double Ck; /* average cost as given by the rule:
	 Qk = Qdecay*Qk + 1, Ck += (fabs (f) - Ck)/Qk */
	double wolfe_hi; /* upper bound for slope in Wolfe test */
	double wolfe_lo; /* lower bound for slope in Wolfe test */
	double awolfe_hi; /* upper bound for slope, approximate Wolfe test */
	int AWolfe; /* F (use Wolfe line search)
	 T (use approximate Wolfe line search)
	 do not change user's AWolfe, this value can be
	 changed based on AWolfeFac */
	int Wolfe; /* T (means code reached the Wolfe part of cg_line */
	double rho; /* either Parm->rho or Parm->nan_rho */
	double alphaold; /* previous value for stepsize alpha */
	double *x; /* current iterate */
	double *xtemp; /* x + alpha*d */
	double *d; /* current search direction */
	double *g; /* gradient at x */
	double *gtemp; /* gradient at x + alpha*d */
	cg_parameter *Parm; /* user parameters */
} cg_com;

/**
 * This class represents an optimization problem. The user should inherit from this class
 * and implement the necessary methods.
 */
class CGDProblem {
public:

	/**
	 * Initializes x.
	 * @param x The variables.
	 * @param n The size.
	 */
	virtual void init(double * x, INT n) = 0;

	/**
	 * Objective function.
	 * @param x The variables.
	 * @param n The size.
	 * @return The value of the objective.
	 */
	virtual double f(double * x, INT n) = 0;

	/**
	 * Gradient.
	 * @param grad_f The gradient (outut).
	 * @param x The variables.
	 * @param n The size.
	 */
	virtual void grad(double * grad_f, double * x, INT n) = 0;

	/**
	 * Compute objective and gradient at the same time. This is a default
	 * implementation that should be overloaded if needed.
	 * @param grad_f The gradient (outut).
	 * @param x The variables.
	 * @param n The size.
	 * @return The value of the objective.
	 */
	virtual double fgrad(double * grad_f, double * x, INT n) {
		grad(grad_f, x, n);
		return f(x, n);
	}

	/**
	 * Finalize the solution.
	 * @param x The variables.
	 * @param n The size.
	 */
	virtual void finalize(double * x, INT n) = 0;

	/**
	 * Destructor.
	 */
	virtual ~CGDProblem() = default;
};

/**
 * The CG algorithm.
 */
class CGDescent {
protected:
	CGDProblem* pb; /**< The optimization problem. */
	double one[1];
	double zero[1];
	BLAS_INT blas_one[1];

	int cg_Wolfe(double alpha, /* stepsize */
	double f, /* function value associated with stepsize alpha */
	double dphi, /* derivative value associated with stepsize alpha */
	cg_com *Com /* cg com */
	);

	int cg_tol(double gnorm, /* gradient sup-norm */
	cg_com *Com /* cg com */
	);

	int cg_line(cg_com *Com /* cg com structure */
	);

	int cg_contract(double *A, /* left side of bracketing interval */
	double *fA, /* function value at a */
	double *dA, /* derivative at a */
	double *B, /* right side of bracketing interval */
	double *fB, /* function value at b */
	double *dB, /* derivative at b */
	cg_com *Com /* cg com structure */
	);

	int cg_evaluate(const char *what, /* fg = evaluate func and grad, g = grad only,f = func only*/
	const char *nan, /* y means check function/derivative values for nan */
	cg_com *Com);

	double cg_cubic(double a, double fa, /* function value at a */
	double da, /* derivative at a */
	double b, double fb, /* function value at b */
	double db /* derivative at b */
	);

	void cg_matvec(double *y, /* product vector */
	double *A, /* dense matrix */
	double *x, /* input vector */
	int n, /* number of columns of A */
	INT m, /* number of rows of A */
	int w /* T => y = A*x, F => y = A'*x */
	);

	void cg_trisolve(double *x, /* right side on input, solution on output */
	double *R, /* dense matrix */
	int m, /* leading dimension of R */
	int n, /* dimension of triangular system */
	int w /* T => Rx = y, F => R'x = y */
	);

	double cg_inf(double *x, /* vector */
	INT n /* length of vector */
	);

	void cg_scale0(double *y, /* output vector */
	double *x, /* input vector */
	double s, /* scalar */
	int n /* length of vector */
	);

	void cg_scale(double *y, /* output vector */
	double *x, /* input vector */
	double s, /* scalar */
	INT n /* length of vector */
	);

	void cg_daxpy0(double *x, /* input and output vector */
	double *d, /* direction */
	double alpha, /* stepsize */
	int n /* length of the vectors */
	);

	void cg_daxpy(double *x, /* input and output vector */
	double *d, /* direction */
	double alpha, /* stepsize */
	INT n /* length of the vectors */
	);

	double cg_dot0(double *x, /* first vector */
	double *y, /* second vector */
	int n /* length of vectors */
	);

	double cg_dot(double *x, /* first vector */
	double *y, /* second vector */
	INT n /* length of vectors */
	);

	void cg_copy0(double *y, /* output of copy */
	double *x, /* input of copy */
	int n /* length of vectors */
	);

	void cg_copy(double *y, /* output of copy */
	double *x, /* input of copy */
	INT n /* length of vectors */
	);

	void cg_step(double *xtemp, /*output vector */
	double *x, /* initial vector */
	double *d, /* search direction */
	double alpha, /* stepsize */
	INT n /* length of the vectors */
	);

	void cg_init(double *x, /* input and output vector */
	double s, /* scalar */
	INT n /* length of vector */
	);

	double cg_update_2(double *gold, /* old g */
	double *gnew, /* new g */
	double *d, /* d */
	INT n /* length of vectors */
	);

	double cg_update_inf(double *gold, /* old g */
	double *gnew, /* new g */
	double *d, /* d */
	INT n /* length of vectors */
	);

	double cg_update_ykyk(double *gold, /* old g */
	double *gnew, /* new g */
	double *Ykyk, double *Ykgk, INT n /* length of vectors */
	);

	double cg_update_inf2(double *gold, /* old g */
	double *gnew, /* new g */
	double *d, /* d */
	double *gnorm2, /* 2-norm of g */
	INT n /* length of vectors */
	);

	double cg_update_d(double *d, double *g, double beta, double *gnorm2, /* 2-norm of g */
	INT n /* length of vectors */
	);

	void cg_Yk(double *y, /*output vector */
	double *gold, /* initial vector */
	double *gnew, /* search direction */
	double *yty, /* y'y */
	INT n /* length of the vectors */
	);

	void cg_printParms(cg_parameter *Parm);

	/**
	 * Set default parameter values.
	 */
	void cg_default(cg_parameter *Parm);

public:
	/**
	 * Constructor.
	 * @param pb The optimization problem to be solved.
	 */
	CGDescent(CGDProblem* pb);

	/**
	 * @param x input: starting guess, output: the solution.
	 * @param n problem dimension.
	 * @param Stat structure with statistics (can be NULL)
	 * @param UParm user parameters, NULL = use default parameters
	 * @param grad_tol StopRule = 1: |g|_infty <= max (grad_tol, StopFac*initial |g|_infty) [default]. StopRule = 0: |g|_infty <= grad_tol(1+|f|)
	 * @param Work NULL => let code allocate memory not NULL => use array Work for required memory. The amount of memory needed depends on the value
	 * of the parameter memory in the Parm structure.  memory > 0 => need (mem+6)*n + (3*mem+9)*mem + 5 where mem = MIN(memory, n) memory = 0 => need 4*n
	 * @return status of solution process:
	 0 (convergence tolerance satisfied)
	 1 (change in func <= feps*|f|)
	 2 (total number of iterations exceeded maxit)
	 3 (slope always negative in line search)
	 4 (number of line search iterations exceeds nline)
	 5 (search direction not a descent direction)
	 6 (excessive updating of eps)
	 7 (Wolfe conditions never satisfied)
	 8 (debugger is on and the function value increases)
	 9 (no cost or gradient improvement in
	 2n + Parm->nslow iterations)
	 10 (out of memory)
	 11 (function nan or +-INF and could not be repaired)
	 12 (invalid choice for memory parameter)
	 * @param init If true the algorithm will request initialization from the problem, otherwise it will take the values initially in x
	 * as the starting point.
	 */
	int cg_descent(double *x, INT n, cg_stats *Stat, cg_parameter *UParm,
			double grad_tol, double *Work, bool init = true);

	/**
	 * Destructor.
	 */
	virtual ~CGDescent() = default;
};
} // namespace
#endif /* CGDESCENT_HPP_ */
