/*
 * cgblas.hpp
 *
 *  Created on: Aug 8, 2016
 *      Author: Said Kerrache
 */

#ifndef CGBLAS_HPP_
#define CGBLAS_HPP_

namespace CG {
/* If the BLAS are not installed, then the following definitions
 can be ignored. If the BLAS are available, then to use them,
 comment out the the next statement (#define NOBLAS) and make
 any needed adjustments to BLAS_UNDERSCORE and the START parameters.
 cg_descent already does loop unrolling, so there is likely no
 benefit from using unrolled BLAS. There could be a benefit from
 using threaded BLAS if the problems is really big. However,
 performing low dimensional operations with threaded BLAS can be
 less efficient than the cg_descent unrolled loops. Hence,
 START parameters should be specified to determine when to start
 using the BLAS. */

#define NOBLAS

/* if BLAS are used, specify the integer precision */
typedef long int BLAS_INT;

/* if BLAS are used, comment out the next statement if no
 *    underscore in the subroutine names are needed */
/* #define BLAS_UNDERSCORE */

/* only use ddot when the vector size >= DDOT_START */
const int DDOT_START = 100;

/* only use dcopy when the vector size >= DCOPY_START */
const int DCOPY_START = 100;

/* only use ddot when the vector size >= DAXPY_START */
const int DAXPY_START = 6000;

/* only use dscal when the vector size >= DSCAL_START */
const int DSCAL_START = 6000;

/* only use idamax when the vector size >= IDAMAX_START */
const int IDAMAX_START = 25;

/* only use matrix BLAS for transpose multiplication when number of
 elements in matrix >= MATVEC_START */
const int MATVEC_START = 8000;

#ifdef BLAS_UNDERSCORE

#define CG_DGEMV dgemv_
#define CG_DTRSV dtrsv_
#define CG_DAXPY daxpy_
#define CG_DDOT ddot_
#define CG_DSCAL dscal_
#define CG_DCOPY dcopy_
#define CG_IDAMAX idamax_

#else

#define CG_DGEMV cblas_dgemv
#define CG_DTRSV cblas_dtrsv
#define CG_DAXPY cblas_daxpy
#define CG_DDOT cblas_ddot
#define CG_DSCAL cblas_dscal
#define CG_DCOPY cblas_dcopy
#define CG_IDAMAX cblas_idamax

#endif

void CG_DGEMV(char *trans, BLAS_INT *m, BLAS_INT *n, double *alpha, double *A,
		BLAS_INT *lda, double *X, BLAS_INT *incx, double *beta, double *Y,
		BLAS_INT *incy);

void CG_DTRSV(char *uplo, char *trans, char *diag, BLAS_INT *n, double *A,
		BLAS_INT *lda, double *X, BLAS_INT *incx);

void CG_DAXPY(BLAS_INT *n, double *DA, double *DX, BLAS_INT *incx, double *DY,
		BLAS_INT *incy);

double CG_DDOT(BLAS_INT *n, double *DX, BLAS_INT *incx, double *DY,
		BLAS_INT *incy);

void CG_DSCAL(BLAS_INT *n, double *DA, double *DX, BLAS_INT *incx);

void CG_DCOPY(BLAS_INT *n, double *DX, BLAS_INT *incx, double *DY,
		BLAS_INT *incy);

BLAS_INT CG_IDAMAX(BLAS_INT *n, double *DX, BLAS_INT *incx);
} //namespace
#endif /* INCLUDE_CGBLAS_HPP_ */
