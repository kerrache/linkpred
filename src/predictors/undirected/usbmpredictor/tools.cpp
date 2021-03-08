/*
 tools.c
 $LastChangedDate: 2011-03-16 13:47:12 +0100 (Wed, 16 Mar 2011) $
 $Revision: 241 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>

#include "linkpred/predictors/undirected/usbmpredictor/tools.h"

/*
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 Vector and matrix memory management
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 */

/*
 ---------------------------------------------------------------------
 Allocation of a matrix of integers
 ---------------------------------------------------------------------
 */
int **
allocate_i_mat(int nrows, int ncolumns) {
	int **array;
	int i;

	array = (int**) malloc(nrows * sizeof(int *));
	if (array == NULL) {
		fprintf(stderr, "out of memory\n");
		return NULL;
	}

	for (i = 0; i < nrows; i++) {
		array[i] = (int*) malloc(ncolumns * sizeof(int));
		if (array[i] == NULL) {
			fprintf(stderr, "out of memory\n");
			return NULL;
		}
	}

	return array;
}

/*
 ---------------------------------------------------------------------
 Allocation of a vector of integers
 ---------------------------------------------------------------------
 */
int *
allocate_i_vec(int nelem) {
	int *array;

	array = (int*) malloc(nelem * sizeof(int));
	if (array == NULL) {
		fprintf(stderr, "out of memory\n");
		return NULL;
	}

	return array;
}

/*
 ---------------------------------------------------------------------
 Allocation of a vector of doubles
 ---------------------------------------------------------------------
 */
double *
allocate_d_vec(int nelem) {
	double *array;

	array = (double*) malloc(nelem * sizeof(double));
	if (array == NULL) {
		fprintf(stderr, "out of memory\n");
		return NULL;
	}

	return array;
}

/*
 ---------------------------------------------------------------------
 Allocation of a matrix of doubles
 ---------------------------------------------------------------------
 */
double **
allocate_d_mat(int nrows, int ncolumns) {
	double **array;
	int i;

	array = (double**) malloc(nrows * sizeof(double *));
	if (array == NULL) {
		fprintf(stderr, "out of memory\n");
		return NULL;
	}

	for (i = 0; i < nrows; i++) {
		array[i] = (double*) malloc(ncolumns * sizeof(double));
		if (array[i] == NULL) {
			fprintf(stderr, "out of memory\n");
			return NULL;
		}
	}

	return array;
}

/*
 ---------------------------------------------------------------------
 Free a matrix of integers
 ---------------------------------------------------------------------
 */
void free_i_mat(int **data, int nrows) {
	int i;

	for (i = 0; i < nrows; i++) {
		free(data[i]);
	}

	free(data);
}

/*
 ---------------------------------------------------------------------
 Free a vector of integers
 ---------------------------------------------------------------------
 */
void free_i_vec(int *data) {
	free(data);
}

/*
 ---------------------------------------------------------------------
 Free a vector of doubles
 ---------------------------------------------------------------------
 */
void free_d_vec(double *data) {
	free(data);
}

/*
 ---------------------------------------------------------------------
 Free a matrix of doubles
 ---------------------------------------------------------------------
 */
void free_d_mat(double **data, int nrows) {
	int i;

	for (i = 0; i < nrows; i++) {
		free(data[i]);
	}

	free(data);
}

/*
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 Random number generation and distributions
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 */

/*
 ---------------------------------------------------------------------
 Returns a number distributed according to the geometric distribution.
 ---------------------------------------------------------------------
 */
int geometric_dist_val(double p, gsl_rng *gen) {
	int val = 0;

	while (gsl_rng_uniform(gen) > p)
		val++;

	return val;
}

/*
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 Statistics
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 */

/*
 ---------------------------------------------------------------------
 Get the mean of the values in an array (of size N).
 ---------------------------------------------------------------------
 */
double mean(double *data, int N) {
	int i;
	double m = 0.0;

	for (i = 0; i < N; i++)
		m += data[i];

	return m / (double) N;
}

/*
 ---------------------------------------------------------------------
 Get the standard deviation of the values in an array (of size N).
 ---------------------------------------------------------------------
 */
double stddev(double *data, int N) {
	int i;
	double m = 0.0, m2 = 0.0;
	double s;

	for (i = 0; i < N; i++) {
		m += data[i];
		m2 += data[i] * data[i];
	}

	m /= (double) (N);
	m2 /= (double) (N);
	if (m2 - m * m > 0.0)
		s = sqrt(m2 - m * m);
	else
		s = 0.0;

	return s;
}

/*
 ---------------------------------------------------------------------
 Get the largest of the values in an array (of size N).
 ---------------------------------------------------------------------
 */
double max(double *data, int N) {
	int i;
	double m = data[0];

	for (i = 1; i < N; i++)
		if (data[i] > m)
			m = data[i];

	return m;
}

/*
 ---------------------------------------------------------------------
 Get the smallest of the values in an array (of size N).
 ---------------------------------------------------------------------
 */
double min(double *data, int N) {
	int i;
	double m = data[0];

	for (i = 1; i < N; i++)
		if (data[i] < m)
			m = data[i];

	return m;
}

/*
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 Other
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 */
/*
 ---------------------------------------------------------------------
 Factorial
 ---------------------------------------------------------------------
 */
long int fact(long int a) {
	if (a == 0)
		return 1;
	else
		return a * fact(a - 1);
}

/*
 ---------------------------------------------------------------------
 Binomial coefficient
 ---------------------------------------------------------------------
 */
long double Choose(int n, int k) {
	int i;
	double accum = 1.0;

	if (k > n)
		return 0.0;

	if (k > n / 2)
		k = n - k; // faster

	for (i = 1; i <= k; i++)
		accum *= (double) (n - k + i) / (double) i;

	return floor(accum + 0.5); // avoid rounding error
}

/*
 ---------------------------------------------------------------------
 Logarithm of the binomial coefficient
 ---------------------------------------------------------------------
 */
double LogChoose(int a, int b) {
	return (double) logl(Choose(a, b));
}

/*
 ---------------------------------------------------------------------
 Initialize a matrix to be used by FastLogChoose
 ---------------------------------------------------------------------
 */
double **
InitializeFastLogChoose(int LogChooseListSize) {
	double **LogChooseList;
	int i, j;

	LogChooseList = allocate_d_mat(LogChooseListSize, LogChooseListSize);
	for (i = 0; i < LogChooseListSize; i++)
		for (j = 0; j < LogChooseListSize; j++)
			LogChooseList[i][j] = -1.0;

	return LogChooseList;
}

/*
 ---------------------------------------------------------------------
 Free a used by FastLogChoose
 ---------------------------------------------------------------------
 */
void FreeFastLogChoose(double **LogChooseList, int LogChooseListSize) {
	free_d_mat(LogChooseList, LogChooseListSize);
	return;
}

/*
 ---------------------------------------------------------------------
 Fast log of the binomial coefficient: checks, first, in a matrix to
 see if the coefficient has been calculated before. At the begining,
 the matrix MUST BE initialized to <0 values (RECOMMENDED: Use
 InitializeFastLogChoose for that purpose).
 ---------------------------------------------------------------------
 */
double FastLogChoose(int r, int l, double **LogChooseList,
		int LogChooseListSize) {
	/*   fprintf(stderr, "Using fast logchoose\n"); */
	if (r < LogChooseListSize) {
		if (LogChooseList[r][l] < 0.) {
			LogChooseList[r][l] = LogChoose(r, l);
		}
		return LogChooseList[r][l];
	} else {
		return LogChoose(r, l);
	}
}

/*
 ---------------------------------------------------------------------
 Initialize a matrix to be used by FastLog
 ---------------------------------------------------------------------
 */
double *
InitializeFastLog(int LogListSize) {
	double *LogList;
	int i;

	LogList = allocate_d_vec(LogListSize);
	for (i = 0; i < LogListSize; i++)
		LogList[i] = log((double) i);

	return LogList;
}

/*
 ---------------------------------------------------------------------
 Free a used by FastLog
 ---------------------------------------------------------------------
 */
void FreeFastLog(double *LogList) {
	free_d_vec(LogList);
	return;
}

/*
 ---------------------------------------------------------------------
 Fast log: if r is small enough, it returns the result from
 previously tabulated values, otherwise, it calculates the
 value. CAUTION!!!! At the begining, the LogList vector MUST BE
 initialized to log(i) values (RECOMMENDED: Use InitializeFastLog for
 that purpose).
 ---------------------------------------------------------------------
 */
double FastLog(int r, double *LogList, int LogListSize) {
	if (r < LogListSize)
		return LogList[r];
	else
		return log((double) r);
}

/*
 ---------------------------------------------------------------------
 Initialize a vector to be used by FastLogFact
 ---------------------------------------------------------------------
 */
double *
InitializeFastLogFact(int LogFactListSize) {
	double *LogFactList;
	int i;

	LogFactList = allocate_d_vec(LogFactListSize);
	for (i = 0; i < LogFactListSize; i++)
		LogFactList[i] = gsl_sf_lnfact(i);

	return LogFactList;
}

/*
 ---------------------------------------------------------------------
 Free a vector used by FastLogFact
 ---------------------------------------------------------------------
 */
void FreeFastLogFact(double *LogFactList) {
	free_d_vec(LogFactList);
	return;
}

/*
 ---------------------------------------------------------------------
 Fast log of the factorial: if r is small enough, it returns the
 result from previously tabulated values, otherwise, it calculates
 the value. CAUTION!!!! At the begining, the LogFactList vector MUST
 BE initialized to log(i!) values (RECOMMENDED: Use
 InitializeFastLogFact for that purpose).
 ---------------------------------------------------------------------
 */
double FastLogFact(int r, double *LogFactList, int LogFactListSize) {
	if (r < LogFactListSize)
		return LogFactList[r];
	else
		return gsl_sf_lnfact(r);
}

/*
 -----------------------------------------------------------------------------
 -----------------------------------------------------------------------------
 File operations
 -----------------------------------------------------------------------------
 -----------------------------------------------------------------------------
 */
/*
 -----------------------------------------------------------------------------
 Count the lines in a file
 -----------------------------------------------------------------------------
 */
int CountLinesInFile(char *inFileName) {
	char line[1000];
	int n = 0;
	FILE *inFile;

	inFile = fopen(inFileName, "r");
	while (fgets(line, 1000, inFile))
		n++;
	fclose(inFile);

	return n;
}
