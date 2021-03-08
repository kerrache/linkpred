/*
 * This file is part of LinkPred.
 *
 * LinkPred: A high performance library for link prediction in complex networks.
 * Copyright (C) 2017  by Said Kerrache.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "linkpred/numerical/linear/gfmatrix.hpp"
#include "linkpred/utils/log.hpp"
#include <limits>
#include <cmath>
#include <cassert>
#include <stdexcept>

//#define LINKPRED_WITH_MKL
#ifdef LINKPRED_WITH_MKL
#include "linkpred/mkl.h"
#endif

namespace LinkPred {


GFMatrix::GFMatrix(int m, int n, bool initZero) {
	logger(logDebug, "GFMatrix constructor...")
	this->m = m;
	this->n = n;
#ifdef LINKPRED_WITH_MKL
	values = (double *) mkl_malloc(m * n * sizeof(double), 64);
#else
	values = new double[m * n];
#endif
	if (initZero) {
#ifdef LINKPRED_WITH_MKL
		double val = 0;
		cblas_dcopy(m * n, &val, 0, values, 1);
#else
		for (int k = 0; k < m * n; k++) {
			values[k] = 0;
		}
#endif
	}
	logger(logDebug, "Done")
}

GFMatrix::GFMatrix(GFMatrix const & that) {
	logger(logDebug, "GFMatrix copy constructor...")
	this->m = that.m;
	this->n = that.n;
#ifdef LINKPRED_WITH_MKL
	values = (double *) mkl_malloc(m * n * sizeof(double), 64);
	cblas_dcopy(m * n, that.values, 1, values, 1);
#else
	values = new double[m * n];
	for (int k = 0; k < m * n; k++) {
		this->values[k] = that.values[k];
	}
#endif
	logger(logDebug, "Done")
}

GFMatrix & GFMatrix::operator =(GFMatrix const & that) {
	logger(logDebug, "GFMatrix copy assignment operator...")
	if (this == &that) {
		return *this;
	}

	if (((m != that.m) || (n != that.n)) && (values != nullptr)) {
#ifdef LINKPRED_WITH_MKL
		mkl_free(values);
		values = (double *) mkl_malloc(m * n * sizeof(double), 64);
#else
		delete[] values;
		values = new double[m * n];
#endif
	}
	m = that.m;
	n = that.n;

#ifdef LINKPRED_WITH_MKL
	cblas_dcopy(m * n, that.values, 1, values, 1);
#else
	double* tv = that.values;
	for (int k = 0; k < m * n; k++) {
		values[k] = tv[k];
	}
#endif
	logger(logDebug, "Done")
	return *this;
}

GFMatrix::GFMatrix(GFMatrix && that) {
	logger(logDebug, "GFMatrix move constructor...")
	m = that.m;
	n = that.n;
	that.m = 0;
	that.n = 0;
	values = that.values;
	that.values = nullptr;
	logger(logDebug, "Done")
}

GFMatrix& GFMatrix::operator =(GFMatrix && that) {
	logger(logDebug, "GFMatrix move assignment operator...")
	if (this == &that) {
		return *this;
	}
	if (values != nullptr) {
#ifdef LINKPRED_WITH_MKL
		mkl_free(values);
#else
		delete[] values;
#endif
	}
	m = that.m;
	n = that.n;
	that.m = 0;
	that.n = 0;
	values = that.values;
	that.values = nullptr;
	logger(logDebug, "Done")
	return *this;
}

bool GFMatrix::operator==(const GFMatrix &that) const {
	logger(logDebug, "GFMatrix ==...")
	if ((this->m != that.m) || (this->n != that.n)) {
		return false;
	}
	for (int k = 0; k < m * n; k++) {
		if (std::abs(this->values[k] - that.values[k])
				>= std::numeric_limits<double>::epsilon()) {
			logger(logDebug, "Done")
			return false;
		}
	}
	logger(logDebug, "Done")
	return true;
}

bool GFMatrix::operator!=(const GFMatrix &that) const {
	return !(*this == that);
}

Vec GFMatrix::getRow(int i) const {
	logger(logDebug, "Getting row " << i << " ...")
	Vec vec(n);

#ifdef LINKPRED_WITH_MKL
	double *vv = vec.values;
	cblas_dcopy(n, &(values[i]), m, vv, 1);
#else
	for (int j = 0; j < n; j++) {
		vec[j] = values[i + j * m];
	}
#endif
	logger(logDebug, "Done")

	return vec;
}

Vec GFMatrix::getCol(int j) const {
	logger(logDebug, "Getting column " << j << " ...")
	Vec vec(m);

#ifdef LINKPRED_WITH_MKL
	double *vv = vec.values;
	cblas_dcopy(m, &(values[j * m]), 1, vv, 1);
#else
	for (int i = 0; i < m; i++) {
		vec[i] = values[i + j * m];
	}
#endif
	logger(logDebug, "Done")
	return vec;
}

void GFMatrix::setRow(int i, Vec const & vec) {
	logger(logDebug, "Setting row " << i << " ...")

#ifdef LINKPRED_WITH_MKL
	double *vv = vec.values;
	cblas_dcopy(n, vv, 1, &(values[i]), m);
#else
	for (int j = 0; j < n; j++) {
		values[i + j * m] = vec[j];
	}
#endif
	logger(logDebug, "Done")
}

void GFMatrix::setCol(int j, Vec const & vec) {
	logger(logDebug, "Setting column " << j << " ...")

#ifdef LINKPRED_WITH_MKL
	double *vv = vec.values;
	cblas_dcopy(m, vv, 1, &(values[j * m]), 1);
#else
	for (int i = 0; i < m; i++) {
		values[i + j * m] = vec[i];
	}
#endif
	logger(logDebug, "Done")
}

GFMatrix GFMatrix::getRows(std::vector<int> const & rowInd) const {
	logger(logDebug, "Getting rows...")
	GFMatrix res(rowInd.size(), n);

#ifdef LINKPRED_WITH_MKL
	double *resv = res.values;
	int nn = rowInd.size();
	int k = 0;
	for (auto it = rowInd.begin(); it != rowInd.end(); ++it, k++) {
		int i = *it;
		cblas_dcopy(n, &(values[i]), m, &(resv[k]), nn);
	}
#else
	double *resv = res.values;
	int nn = rowInd.size();
	int k = 0;
	for (auto it = rowInd.begin(); it != rowInd.end(); ++it, k++) {
		int i = *it;
		for (int j = 0; j < n; j++) {
			resv[k + j * nn] = values[i + j * m];
		}
	}
#endif
	logger(logDebug, "Done")
	return res;
}

GFMatrix GFMatrix::getCols(std::vector<int> const & colInd) const {
	logger(logDebug, "Getting columns...")
	GFMatrix res(m, colInd.size());

#ifdef LINKPRED_WITH_MKL
	double *resv = res.values;
	int k = 0;
	for (auto it = colInd.begin(); it != colInd.end(); ++it, k++) {
		int j = *it;
		cblas_dcopy(m, &(values[j * m]), 1, &(resv[k * m]), 1);
	}
#else
	double *resv = res.values;
	int k = 0;
	for (auto it = colInd.begin(); it != colInd.end(); ++it) {
		int j = *it;
		for (int i = 0; i < m; i++, k++) {
			resv[k] = values[i + j * m];
		}
	}
#endif
	logger(logDebug, "Done")
	return res;
}
Vec operator *(const GFMatrix &mat, const Vec &vec) {

	logger(logDebug, "Matrix * Vec...")

	if (mat.n != vec.n) {
		throw std::invalid_argument("Incompatible matrix dimensions");
	}

	Vec res(mat.m);

#ifdef LINKPRED_WITH_MKL
	cblas_dgemv(CblasColMajor, CblasNoTrans, mat.m, mat.n, 1, mat.values, mat.m,
			vec.values, 1, 0, res.values, 1);
#else
	for (int i = 0; i < mat.m; i++) {
		double s = 0;
		for (int j = 0; j < mat.n; j++) {
			logger(logDebug1, i << " " << j << " " << mat.get(i, j))
			s += mat.values[i + j * mat.m] * vec[j];
		}
		res[i] = s;
	}
#endif
	logger(logDebug, "Done")
	return res;
}

Vec GFMatrix::sumCols() const {
	logger(logDebug, "Summing matrix columns...")
	Vec sum(n);
	int k = 0;
	for (int j = 0; j < n; j++) {
		double s = 0;
		for (int i = 0; i < m; i++, k++) {
			s += values[k];
		}
		sum[j] = s;
	}
	logger(logDebug, "Done")
	return sum;
}

Vec GFMatrix::sumRows() const {
	logger(logDebug, "Summing matrix rows...")
	Vec sum(m);
	for (int i = 0; i < m; i++) {
		sum[i] = 0;
	}
	int k = 0;
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++, k++) {
			sum[i] += values[k];
		}
	}
	logger(logDebug, "Done")
	return sum;
}

double GFMatrix::sum() {
	logger(logDebug, "Summing all matrix entries...")
	double s = 0;
	for (int k = 0; k < m * n; k++) {
		s += values[k];
	}
	logger(logDebug, "Done")
	return s;
}

void GFMatrix::print() const {
	logger(logDebug, "Printing matrix...")
	std::cout << std::endl;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << get(i, j) << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	logger(logDebug, "Done")
}

void GFMatrix::print(std::string name) const {
	logger(logDebug, "Printing matrix...")
	std::cout << name << ": " << std::endl;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << get(i, j) << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	logger(logDebug, "Done")
}

Vec GFMatrix::diag() const {
	logger(logDebug, "Extracting matrix diagonal...")
	if (n != m) {
		throw std::invalid_argument(
				"Matrix is not symmetric. Cannot extract diagonal");
	}
	Vec diag(n);
	for (int i = 0; i < n; i++) {
		diag[i] = get(i, i);
	}
	logger(logDebug, "Done")
	return diag;
}

void GFMatrix::setDiag(double v) {
	logger(logDebug, "Setting matrix diagonal...")
	if (n != m) {
		throw std::invalid_argument(
				"Matrix is not symmetric. Cannot extract diagonal. Abort");
	}
	for (int i = 0; i < n; i++) {
		set(i, i, v);
	}
	logger(logDebug, "Done")
}

void GFMatrix::removeNaN() {
	for (int i = 0; i < m * n; i++) {
		if (values[i] != values[i]) { // Testing for NaN
			values[i] = 0;
		}
	}
}

GFMatrix::~GFMatrix() {
	logger(logDebug, "GFMatrix " << name << " destructor...")
	if (values != nullptr) {
#ifdef LINKPRED_WITH_MKL
		mkl_free(values);
#else
		delete[] values;
#endif
	}
	logger(logDebug, "Done")
}

void GFMatrix::matMultFF(double const *values1, double const *values2,
		double * values3, int m1, int n1, int m2, int n2) {
	if (n1 != m2) {
		throw std::invalid_argument("Incompatible matrix dimensions");
	}

	for (int j = 0; j < n2; j++) {
		for (int i = 0; i < m1; i++) {
			double sum = 0;
			for (int l = 0; l < n1; l++) {
				sum += values1[i + l * m1] * values2[l + j * m2];
			}
			values3[i + j * m1] = sum;
		}
	}
}

void GFMatrix::matMultFT(double const *values1, double const *values2,
		double * values3, int m1, int n1, int m2, int n2) {
	if (n1 != n2) {
		throw std::invalid_argument("Incompatible matrix dimensions");
	}

	for (int j = 0; j < m2; j++) {
		for (int i = 0; i < m1; i++) {
			double sum = 0;
			for (int l = 0; l < n1; l++) {
				sum += values1[i + l * m1] * values2[j + l * m2];
			}
			values3[i + j * m1] = sum;
		}
	}
}

void GFMatrix::matMultTF(double const *values1, double const *values2,
		double * values3, int m1, int n1, int m2, int n2) {
	if (m1 != m2) {
		throw std::invalid_argument("Incompatible matrix dimensions");
	}

	for (int j = 0; j < n2; j++) {
		for (int i = 0; i < n1; i++) {
			double sum = 0;
			for (int l = 0; l < m1; l++) {
				sum += values1[l + i * m1] * values2[l + j * m2];
			}
			values3[i + j * n1] = sum;
		}
	}
}

void GFMatrix::matMultTT(double const *values1, double const *values2,
		double * values3, int m1, int n1, int m2, int n2) {
	if (m1 != n2) {
		throw std::invalid_argument("Incompatible matrix dimensions");
	}

	for (int j = 0; j < m2; j++) {
		for (int i = 0; i < n1; i++) {
			double sum = 0;
			for (int l = 0; l < m1; l++) {
				sum += values1[l + i * m1] * values2[j + l * m2];
			}
			values3[i + j * n1] = sum;
		}
	}
}

GFMatrix GFMatrix::mult(GFMatrix const & mat1, GFMatrix const & mat2,
		bool trans1, bool trans2) {
	logger(logDebug, "Matrix * Matrix...")

#ifdef LINKPRED_WITH_MKL
	CBLAS_TRANSPOSE transa = (trans1) ? CblasTrans : CblasNoTrans;
	CBLAS_TRANSPOSE transb = (trans2) ? CblasTrans : CblasNoTrans;

	if (trans1) {
		if (trans2) {
			GFMatrix mat3(mat1.n, mat2.m);
			cblas_dgemm(CblasColMajor, transa, transb, mat1.n, mat2.m, mat1.m,
					1, mat1.values, mat1.m, mat2.values, mat2.m, 0, mat3.values,
					mat3.m);
			return mat3;
		} else {
			GFMatrix mat3(mat1.n, mat2.n);
			cblas_dgemm(CblasColMajor, transa, transb, mat1.n, mat2.n, mat1.m,
					1, mat1.values, mat1.m, mat2.values, mat2.m, 0, mat3.values,
					mat3.m);
			return mat3;
		}
	} else {
		if (trans2) {
			GFMatrix mat3(mat1.m, mat2.m);
			cblas_dgemm(CblasColMajor, transa, transb, mat1.m, mat2.m, mat1.n,
					1, mat1.values, mat1.m, mat2.values, mat2.m, 0, mat3.values,
					mat3.m);
			return mat3;

		} else {
			GFMatrix mat3(mat1.m, mat2.n);
			cblas_dgemm(CblasColMajor, transa, transb, mat1.m, mat2.n, mat1.n,
					1, mat1.values, mat1.m, mat2.values, mat2.m, 0, mat3.values,
					mat3.m);
			return mat3;
		}
	}
#else
	if (trans1) {
		if (trans2) {
			GFMatrix mat3(mat1.n, mat2.m);
			matMultTT(mat1.values, mat2.values, mat3.values, mat1.m, mat1.n,
					mat2.m, mat2.n);
			return mat3;
		} else {
			GFMatrix mat3(mat1.n, mat2.n);
			matMultTF(mat1.values, mat2.values, mat3.values, mat1.m, mat1.n,
					mat2.m, mat2.n);
			return mat3;
		}
	} else {
		if (trans2) {
			GFMatrix mat3(mat1.m, mat2.m);
			matMultFT(mat1.values, mat2.values, mat3.values, mat1.m, mat1.n,
					mat2.m, mat2.n);
			return mat3;

		} else {
			GFMatrix mat3(mat1.m, mat2.n);
			matMultFF(mat1.values, mat2.values, mat3.values, mat1.m, mat1.n,
					mat2.m, mat2.n);
			return mat3;
		}
	}
#endif
	logger(logDebug, "Done")
}

GFMatrix GFMatrix::elemMult(GFMatrix const & mat1, GFMatrix const & mat2) {
	logger(logDebug, "Matrix .* Matrix...")
	if ((mat1.n != mat2.n) || (mat1.m != mat2.m)) {
		throw std::invalid_argument("Incompatible matrix dimensions");
	}

	int m = mat1.m;
	int n = mat1.n;
	GFMatrix mat3(m, n);
	double *v1 = mat1.values;
	double *v2 = mat2.values;
	double *v3 = mat3.values;
#ifdef LINKPRED_WITH_MKL
	vdMul(m * n, v1, v2, v3);
#else
	for (int i = 0; i < m * n; i++) {
		v3[i] = v1[i] * v2[i];
	}
#endif
	logger(logDebug, "Done")
	return mat3;
}

GFMatrix GFMatrix::mult(Vec const & v1, Vec const & v2) {
	logger(logDebug, "Vec' * Vec...")
	int m = v1.size();
	int n = v2.size();
	GFMatrix res(m, n);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			res.set(i, j, v1[i] * v2[j]);
		}
	}
	logger(logDebug, "Done")
	return res;
}

GFMatrix operator*(const GFMatrix & mat1, const GFMatrix & mat2) {
	logger(logDebug, "Matrix * Matrix...")
#ifdef LINKPRED_WITH_MKL
	int m = mat1.m;
	int k = mat1.n;
	if (k != mat2.m) {
		throw std::invalid_argument("Incompatible matrix dimensions");
	}
	int n = mat2.n;
	GFMatrix mat3(m, n);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1,
			mat1.values, m, mat2.values, k, 0, mat3.values, m);

#else
	GFMatrix mat3(mat1.m, mat2.n);
	auto values1 = mat1.values;
	auto values2 = mat2.values;
	auto values3 = mat3.values;
	GFMatrix::matMultFF(values1, values2, values3, mat1.m, mat1.n, mat2.m,
			mat2.n);
#endif

	logger(logDebug, "Done")
	return mat3;
}

GFMatrix operator+(const GFMatrix & mat1, const GFMatrix & mat2) {
	logger(logDebug, "Matrix + Matrix...")

	if ((mat1.n != mat2.n) || (mat1.m != mat2.m)) {
		throw std::invalid_argument("Incompatible matrix dimensions");
	}

	int m = mat1.m;
	int n = mat1.n;
	GFMatrix mat3(m, n);
	double *v1 = mat1.values;
	double *v2 = mat2.values;
	double *v3 = mat3.values;
#ifdef LINKPRED_WITH_MKL
	vdAdd(m * n, v1, v2, v3);
#else
	for (int i = 0; i < m * n; i++) {
		v3[i] = v1[i] + v2[i];
	}
#endif
	logger(logDebug, "Done")
	return mat3;
}

GFMatrix operator/(const GFMatrix & mat1, const GFMatrix & mat2) {
	logger(logDebug, "Matrix ./ Matrix...")
	if ((mat1.n != mat2.n) || (mat1.m != mat2.m)) {
		throw std::invalid_argument("Incompatible matrix dimensions");
	}

	int m = mat1.m;
	int n = mat1.n;
	GFMatrix mat3(m, n);
	double *v1 = mat1.values;
	double *v2 = mat2.values;
	double *v3 = mat3.values;
#ifdef LINKPRED_WITH_MKL
	vdDiv(m * n, v1, v2, v3);
#else
	for (int i = 0; i < m * n; i++) {
		v3[i] = v1[i] / v2[i];
	}
#endif
	logger(logDebug, "Done")
	return mat3;
}

GFMatrix operator+(const GFMatrix & mat1, const Vec & mat2) {
	logger(logDebug, "Matrix + Diagonal...")
	if ((mat1.m != mat1.n) || (mat1.m != mat2.n)) {
		throw std::invalid_argument("Incompatible matrix dimensions");
	}

	int n = mat1.n;
	GFMatrix mat3(n, n);
	double *v1 = mat1.values;
	double *v2 = mat2.values;
	double *v3 = mat3.values;
	for (int i = 0; i < n * n; i++) {
		v3[i] = v1[i];
	}
	for (int i = 0; i < n; i++) {
		v3[i + i * n] += v2[i];
	}
	logger(logDebug, "Done")
	return mat3;
}

GFMatrix operator+(const Vec & mat1, const GFMatrix & mat2) {
	return mat2 + mat1;
}

GFMatrix operator-(const GFMatrix & mat1, const Vec & mat2) {
	logger(logDebug, "Matrix - Diagonal...")
	if ((mat1.m != mat1.n) || (mat1.m != mat2.n)) {
		throw std::invalid_argument("Incompatible matrix dimensions");
	}

	int n = mat1.n;
	GFMatrix mat3(n, n);
	double *v1 = mat1.values;
	double *v2 = mat2.values;
	double *v3 = mat3.values;
	for (int i = 0; i < n * n; i++) {
		v3[i] = v1[i];
	}
	for (int i = 0; i < n; i++) {
		v3[i + i * n] -= v2[i];
	}
	logger(logDebug, "Done")
	return mat3;
}

GFMatrix operator-(const Vec & mat1, const GFMatrix & mat2) {
	logger(logDebug, "Diagonal + Matrix...")

	if ((mat2.m != mat2.n) || (mat2.m != mat1.n)) {
		throw std::invalid_argument("Incompatible matrix dimensions");
	}

	int n = mat1.n;
	GFMatrix mat3(n, n);
	double *v1 = mat1.values;
	double *v2 = mat2.values;
	double *v3 = mat3.values;
	for (int i = 0; i < n * n; i++) {
		v3[i] = -v1[i];
	}
	for (int i = 0; i < n; i++) {
		v3[i + i * n] += v2[i];
	}
	logger(logDebug, "Done")
	return mat3;
}

GFMatrix operator*(double a, const GFMatrix & mat) {
	logger(logDebug, "a * Matrix...")
	int m = mat.m;
	int n = mat.n;
	GFMatrix res(m, n);
	double *v1 = mat.values;
	double *v2 = res.values;
	for (int i = 0; i < m * n; i++) {
		v2[i] = a * v1[i];
	}
	logger(logDebug, "Done")
	return res;
}

GFMatrix operator*(const GFMatrix & mat, double a) {
	logger(logDebug, "Matrix * a...")
	int m = mat.m;
	int n = mat.n;
	GFMatrix res(m, n);
	double *v1 = mat.values;
	double *v2 = res.values;
	for (int i = 0; i < m * n; i++) {
		v2[i] = a * v1[i];
	}
	logger(logDebug, "Done")
	return res;
}

GFMatrix operator+(double a, const GFMatrix & mat) {
	logger(logDebug, "a + Matrix...")
	int m = mat.m;
	int n = mat.n;
	GFMatrix res(m, n);
	double *v1 = mat.values;
	double *v2 = res.values;
	for (int i = 0; i < m * n; i++) {
		v2[i] = a + v1[i];
	}
	logger(logDebug, "Done")
	return res;
}

GFMatrix operator+(const GFMatrix & mat, double a) {
	logger(logDebug, "Matrix + a...")
	int m = mat.m;
	int n = mat.n;
	GFMatrix res(m, n);
	double *v1 = mat.values;
	double *v2 = res.values;
	for (int i = 0; i < m * n; i++) {
		v2[i] = v1[i] + a;
	}
	logger(logDebug, "Done")
	return res;
}

GFMatrix operator-(double a, const GFMatrix & mat) {
	logger(logDebug, "a - Matrix...")
	int m = mat.m;
	int n = mat.n;
	GFMatrix res(m, n);
	double *v1 = mat.values;
	double *v2 = res.values;
	for (int i = 0; i < m * n; i++) {
		v2[i] = a - v1[i];
	}
	logger(logDebug, "Done")
	return res;
}

GFMatrix operator-(const GFMatrix & mat, double a) {
	logger(logDebug, "Matrix + a...")
	int m = mat.m;
	int n = mat.n;
	GFMatrix res(m, n);
	double *v1 = mat.values;
	double *v2 = res.values;
	for (int i = 0; i < m * n; i++) {
		v2[i] = v1[i] - a;
	}
	logger(logDebug, "Done")
	return res;
}

GFMatrix operator-(const GFMatrix & mat1, const GFMatrix & mat2) {
	logger(logDebug, "Matrix - Marix...")
	if ((mat1.n != mat2.n) || (mat1.m != mat2.m)) {
		throw std::invalid_argument("Incompatible matrix dimensions");
	}

	int m = mat1.m;
	int n = mat1.n;
	GFMatrix mat3(m, n);
	double *v1 = mat1.values;
	double *v2 = mat2.values;
	double *v3 = mat3.values;
#ifdef LINKPRED_WITH_MKL
	vdSub(m * n, v1, v2, v3);
#else
	for (int i = 0; i < m * n; i++) {
		v3[i] = v1[i] - v2[i];
	}
#endif
	logger(logDebug, "Done")
	return mat3;
}


} /* namespace LinkPred */
