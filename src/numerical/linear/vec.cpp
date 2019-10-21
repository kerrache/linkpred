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

#include "linkpred/numerical/linear/vec.hpp"
#include "linkpred/utils/log.hpp"

#include <fstream>
#include <cmath>
#include <limits>
#include <cfloat>
#include <cstdlib>

#ifdef WITH_MKL
#include "linkpred/mkl.h"
#endif

namespace LinkPred {
// TODO change new to std::malloc in order to use std::realloc
Vec::Vec() {
	n = 0;
	values = nullptr;
	constant = false;
}

Vec::Vec(int _n) {
	n = _n;
#ifdef WITH_MKL
	values = (double *) mkl_malloc(n * sizeof(double), 64);
#else
	values = new double[n];
#endif
	constant = false;
}

Vec::Vec(int n, double * vals) {
	this->n = n;
#ifdef WITH_MKL
	values = (double *) mkl_malloc(n * sizeof(double), 64);
	cblas_dcopy(n, vals, 1, values, 1);
#else
	values = new double[n];
	for (int i = 0; i < n; i++) {
		values[i] = vals[i];
	}
#endif
	constant = false;
}

Vec::Vec(Vec const & v, std::shared_ptr<std::vector<int>> ind) {
	n = ind->size();
#ifdef WITH_MKL
	values = (double *) mkl_malloc(n * sizeof(double), 64);
	cblas_dgthr(n, v.values, values, ind->data());
#else
	values = new double[n];
	for (int i = 0; i < n; i++) {
		values[i] = v.values[(*ind)[i]];
	}
#endif
	constant = false;
}

Vec::Vec(Vec const & v1, Vec const & v2) {
	n = v1.n + v2.n;
#ifdef WITH_MKL
	values = (double *) mkl_malloc(n * sizeof(double), 64);
	cblas_dcopy(v1.n, v1.values, 1, values, 1);
	cblas_dcopy(v2.n, v2.values, 1, &(values[v1.n]), 1);
#else
	values = new double[n];
	double* vv1 = v1.values;
	for (int i = 0; i < v1.n; i++) {
		values[i] = vv1[i];
	}
	double* vv2 = v2.values;
	for (int i = 0; i < v2.n; i++) {
		values[i + v1.n] = vv2[i];
	}
#endif
	constant = false;
}

Vec::Vec(Vec const & that) {
	logger(logDebug, "Vec copy constructor...")
	n = that.n;
#ifdef WITH_MKL
	values = (double *) mkl_malloc(n * sizeof(double), 64);
	cblas_dcopy(n, that.values, 1, values, 1);
#else
	values = new double[n];
	double* tv = that.values;
	for (int i = 0; i < n; i++) {
		values[i] = tv[i];
	}
#endif
	constant = that.constant;
	logger(logDebug, "Done")
}

Vec& Vec::operator =(Vec const & that) {
	logger(logDebug, "Vec copy assignment operator...")

	if (this == &that) {
		logger(logDebug, "Done")
		return *this;
	}

	if (n != that.n) {
		n = that.n;
#ifdef WITH_MKL
		if (values != nullptr) {
			mkl_free(values);
		}
		values = (double *) mkl_malloc(n * sizeof(double), 64);
#else
		if (values != nullptr) {
			delete[] values;
		}
		values = new double[n];
#endif
	}

#ifdef WITH_MKL
	cblas_dcopy(n, that.values, 1, values, 1);
#else
	double* tv = that.values;
	for (int i = 0; i < n; i++) {
		values[i] = tv[i];
	}
#endif
	constant = that.constant;
	logger(logDebug, "Done")
	return *this;
}

Vec::Vec(Vec && that) {
	logger(logDebug, "Vec move constructor...")
	n = that.n;
	that.n = 0;
	values = that.values;
	that.values = nullptr;
	constant = that.constant;
	that.constant = false;
	logger(logDebug, "Done")
}

Vec& Vec::operator =(Vec && that) {
	logger(logDebug, "Vec move assignment operator...")
	if (this == &that) {
		logger(logDebug, "Done")
		return *this;
	}
	if (values != nullptr) {
#ifdef WITH_MKL
		mkl_free(values);
#else
		delete[] values;
#endif
	}
	n = that.n;
	that.n = 0;
	values = that.values;
	that.values = nullptr;
	constant = that.constant;
	that.constant = false;
	logger(logDebug, "Done")
	return *this;
}

bool Vec::operator==(const Vec &that) const {
	logger(logDebug, "Comparing vectors")

	if (n != that.n) {
		logger(logDebug3, "Different dimensions")
		return false;
	}

	for (int k = 0; k < n; k++) {
		if (fabs(values[k] - that.values[k])
				>= std::numeric_limits<double>::epsilon()) {
			logger(logDebug3,
					"Vectors have different values at " + std::to_string(k)
							+ " : " + std::to_string(values[k]) + " - "
							+ std::to_string(that.values[k]))
			return false;
		}
	}

	logger(logDebug3, "Vectors equal")

	return true;
}

bool Vec::operator!=(const Vec &that) const {
	return !(*this == that);
}

void Vec::resize(int newSize) {

	if (newSize > n) {
#ifdef WITH_MKL
		values = (double *) mkl_realloc(values, newSize * sizeof(double));
		if (values == nullptr) {
			logger(logError, "Cannot reallocate memory. Abort.")
			exit(1);
		}
#else
		double * newValues = new double[newSize];
		for (int i = 0; i < n; i++) {
			newValues[i] = values[i];
		}
		delete[] values;
		values = newValues;
#endif
		n = newSize;
	}
}

bool Vec::makeConstant(double tol) {
	double val = values[0];
	if (val != 0) {
		for (int i = 1; i < n; i++) {
			if (std::abs(values[i] / val - 1) > tol) {
				return false;
			}
		}
	} else {
		for (int i = 1; i < n; i++) {
			if (std::abs(values[i]) > tol) {
				return false;
			}
		}
	}
	constant = true;
	return true;
}

void Vec::resize(int newSize, double val) {

	if (newSize > n) {
		double* oldV = values;
#ifdef WITH_MKL
		values = (double *) mkl_malloc(newSize * sizeof(double), 64);
		cblas_dcopy(n, oldV, 1, values, 1);
		mkl_free(oldV);
		cblas_dcopy(newSize - n, &val, 0, &(values[n]), 1); // Set incx = 0 to copy val everywhere
#else
		values = new double[newSize];
		for (int i = 0; i < n; i++) {
			values[i] = oldV[i];
		}
		delete[] oldV;
		for (int i = n; i < newSize; i++) {
			values[i] = val;
		}
#endif
		n = newSize;
	}
}

std::vector<int> Vec::findNz() const {
	std::vector<int> nz;
	for (int i = 0; i < n; i++) {
		if (values[i] != 0) {
			nz.push_back(i);
		}
	}
	return nz;
}

std::vector<int> Vec::find(double val) const {
	std::vector<int> ind;
	for (int i = 0; i < n; i++) {
		if (values[i] == val) {
			ind.push_back(i);
		}
	}
	return ind;
}

void Vec::reset() {
#ifdef WITH_MKL
	double val = 0;
	cblas_dcopy(n, &val, 0, values, 1);
#else
	for (int i = 0; i < n; i++) {
		values[i] = 0;
	}
#endif
}

void Vec::reset(double val) {
#ifdef WITH_MKL
	cblas_dcopy(n, &val, 0, values, 1);
#else
	for (int i = 0; i < n; i++) {
		values[i] = val;
	}
#endif
}

void Vec::read(std::ifstream *in) {

/// read data
	for (int k = 0; k < n; k++) {
		(*in) >> values[k];
	}
}

void Vec::write(std::ofstream *out) const {

// Set precision
	out->precision(std::numeric_limits<double>::digits);

/// write data
	for (int k = 0; k < n; k++) {
		(*out) << values[k] << std::endl;
	}
}

void Vec::read(std::string fileName) {
	std::ifstream in;
	in.open(fileName.c_str(), std::fstream::in);
	if (!in) {
		std::cerr << "Cannot open file: " << fileName.c_str() << std::endl;
		exit(1);
	}
	read(&in);
	in.close();
}

void Vec::write(std::string fileName) const {
	std::ofstream out;
	out.open(fileName.c_str(), std::fstream::out);
	if (!out) {
		std::cerr << "Cannot open file: " << fileName.c_str() << std::endl;
		exit(1);
	}
	write(&out);
	out.close();
}

double Vec::norm() const {
#ifdef WITH_MKL
	return cblas_dnrm2(n, values, 1) / sqrt(n);
#else
	double sum = 0;
	for (int i = 0; i < n; i++) {
		double val = values[i];
		sum += val * val;
	}
	return sqrt(sum / n);
#endif
}

double Vec::mean() const {

	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += values[i];
	}
	return sum / n;
}

double Vec::amin() const {
#ifdef WITH_MKL
	return values[cblas_idamin(n, values, 1)];
#else
	double m = values[0];
	for (int i = 1; i < n; i++) {
		if (std::abs(values[i]) < m) {
			m = values[i];
		}
	}
	return m;
#endif
}

double Vec::amax() const {
#ifdef WITH_MKL
	return values[cblas_idamax(n, values, 1)];
#else
	double m = values[0];
	for (int i = 1; i < n; i++) {
		if (std::abs(values[i]) > m) {
			m = values[i];
		}
	}
	return m;
#endif
}

double Vec::min() const {
	double m = values[0];
	for (int i = 1; i < n; i++) {
		if (values[i] < m) {
			m = values[i];
		}
	}
	return m;
}

double Vec::max() const {
	double m = values[0];
	for (int i = 1; i < n; i++) {
		if (values[i] > m) {
			m = values[i];
		}
	}
	return m;
}

int Vec::ifmin() const {
	double m = values[0];
	int im = 0;
	for (int i = 1; i < n; i++) {
		if (values[i] < m) {
			m = values[i];
			im = i;
		}
	}
	return im;
}

int Vec::ifmax() const {
	double m = values[0];
	int im = 0;
	for (int i = 1; i < n; i++) {
		if (values[i] > m) {
			m = values[i];
			im = i;
		}
	}
	return im;
}

int Vec::ilmin() const {
	double m = values[0];
	int im = 0;
	for (int i = 1; i < n; i++) {
		if (values[i] <= m) {
			m = values[i];
			im = i;
		}
	}
	return im;
}

int Vec::ilmax() const {
	double m = values[0];
	int im = 0;
	for (int i = 1; i < n; i++) {
		if (values[i] >= m) {
			m = values[i];
			im = i;
		}
	}
	return im;
}

void Vec::center() {
	double m = mean();
	for (int i = 0; i < n; i++) {
		values[i] -= m;
	}

}

double Vec::norm(int lim) const {
#ifdef WITH_MKL
	return cblas_dnrm2(lim, values, 1) / sqrt(lim);
#else
	double sum = 0;
	for (int i = 0; i < lim; i++) {
		double val = values[i];
		sum += val * val;
	}
	return sqrt(sum / lim);
#endif

}

double Vec::mean(int lim) const {

	double sum = 0;
	for (int i = 0; i < lim; i++) {
		sum += values[i];
	}
	return sum / lim;
}

double Vec::amin(int lim) const {
#ifdef WITH_MKL
	return values[cblas_idamin(lim, values, 1)];
#else
	double m = values[0];
	for (int i = 1; i < lim; i++) {
		if (std::abs(values[i]) < m) {
			m = values[i];
		}
	}
	return m;
#endif
}

double Vec::amax(int lim) const {
#ifdef WITH_MKL
	return values[cblas_idamax(lim, values, 1)];
#else
	double m = values[0];
	for (int i = 1; i < lim; i++) {
		if (std::abs(values[i]) > m) {
			m = values[i];
		}
	}
	return m;
#endif
}

double Vec::min(int lim) const {
	double m = values[0];
	for (int i = 1; i < lim; i++) {
		if (values[i] < m) {
			m = values[i];
		}
	}
	return m;
}

double Vec::max(int lim) const {
	double m = values[0];
	for (int i = 1; i < lim; i++) {
		if (values[i] > m) {
			m = values[i];
		}
	}
	return m;
}

void Vec::center(int lim) {
	double m = mean(lim);
	for (int i = 0; i < lim; i++) {
		values[i] -= m;
	}
}

void Vec::axpy(double a, Vec const & x) {
	if (x.size() != n) {
		throw std::invalid_argument(
				"Nonconformant arguments: vectors must have the same size.");
	}

#ifdef WITH_MKL
	cblas_daxpy(n, a, x.values, 1, values, 1);
#else
	double*xv = x.values;
	for (int i = 0; i < n; i++) {
		values[i] += a * xv[i];
	}
#endif
}

void Vec::scale(double a) {
#ifdef WITH_MKL
	cblas_dscal(n, a, values, 1);
#else
	for (int i = 0; i < n; i++) {
		values[i] *= a;
	}
#endif
}

void Vec::partScale(int start, int nb, double a) {
#ifdef WITH_MKL
	cblas_dscal(nb, a, &(values[start]), 1);
#else
	for (int i = start; i < start + nb; i++) {
		values[i] *= a;
	}
#endif
}

void Vec::scatter(int start, Vec const & v, int const * ind) {
//	for (int i = 0; i < v.n; i++) {
//		logger(logError, "ind: " << ind[i])
//	}
#ifdef WITH_MKL
	cblas_dsctr(v.n, v.values, ind, &(values[start]));
#else
	double* vv1 = &(values[start]);
	double* vv2 = v.values;
	for (int i = 0; i < v.n; i++) {
		vv1[ind[i]] = vv2[i];
	}
#endif

}

double Vec::sum() const {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += values[i];
	}
	return sum;
}

void Vec::print() const {
	print(std::cout);
}

void Vec::print(std::string header) const {
	print(std::cout, header);
}

void Vec::print(std::ostream & out) const {
	for (int i = 0; i < n; i++) {
		out << values[i] << std::endl;
	}
}

void Vec::print(std::ostream & out, std::string header) const {
	out << header << " :" << std::endl;
	print(out);
}

Vec::~Vec() {
	if (values != nullptr) {
#ifdef WITH_MKL
		mkl_free(values);
#else
		delete[] values;
#endif
	}
}

Vec operator +(Vec const & v1, Vec const & v2) {
	logger(logDebug, "Vec + Vec")
	if (v1.n != v2.n) {
		throw std::invalid_argument(
				"Nonconformant arguments: vectors must have the same size.");
	}
	if (v1.isConstant()) {
		return v1[0] + v2;
	}
	if (v2.isConstant()) {
		return v1 + v2[0];
	}
	int n = v1.n;
	Vec v3(n);
	double *vv1 = v1.values;
	double *vv2 = v2.values;
	double *vv3 = v3.values;
#ifdef WITH_MKL
	vdAdd(n, vv1, vv2, vv3);
#else
	for (int i = 0; i < n; i++) {
		vv3[i] = vv1[i] + vv2[i];
	}
#endif
	logger(logDebug, "Done")
	return v3;
}

Vec operator -(Vec const & v1, Vec const & v2) {
	logger(logDebug, "Vec - Vec")
	if (v1.n != v2.n) {
		throw std::invalid_argument(
				"Nonconformant arguments: vectors must have the same size.");
	}
	if (v1.isConstant()) {
		return v1[0] - v2;
	}
	if (v2.isConstant()) {
		return v1 - v2[0];
	}

	int n = v1.n;
	Vec v3(n);
	double *vv1 = v1.values;
	double *vv2 = v2.values;
	double *vv3 = v3.values;
#ifdef WITH_MKL
	vdSub(n, vv1, vv2, vv3);
#else
	for (int i = 0; i < n; i++) {
		vv3[i] = vv1[i] - vv2[i];
	}
#endif
	logger(logDebug, "Done")
	return v3;
}

Vec operator *(Vec const & v1, Vec const & v2) {
	logger(logDebug, "Vec * Vec")
	if (v1.n != v2.n) {
		throw std::invalid_argument(
				"Nonconformant arguments: vectors must have the same size.");
	}
	if (v1.isConstant()) {
		return v1[0] * v2;
	}
	if (v2.isConstant()) {
		return v1 * v2[0];
	}

	int n = v1.n;
	Vec v3(n);
	double *vv1 = v1.values;
	double *vv2 = v2.values;
	double *vv3 = v3.values;
#ifdef WITH_MKL
	vdMul(n, vv1, vv2, vv3);
#else
	for (int i = 0; i < n; i++) {
		vv3[i] = vv1[i] * vv2[i];
	}
#endif
	logger(logDebug, "Done")
	return v3;
}

Vec operator/(Vec const & v1, Vec const & v2) {
	logger(logDebug, "Vec / Vec")
	if (v1.n != v2.n) {
		throw std::invalid_argument(
				"Nonconformant arguments: vectors must have the same size.");
	}
	if (v1.isConstant()) {
		return v1[0] / v2;
	}
	if (v2.isConstant()) {
		return v1 / v2[0];
	}

	int n = v1.n;
	Vec v3(n);
	double *vv1 = v1.values;
	double *vv2 = v2.values;
	double *vv3 = v3.values;
#ifdef WITH_MKL
	vdDiv(n, vv1, vv2, vv3);
#else
	for (int i = 0; i < n; i++) {
		vv3[i] = vv1[i] / vv2[i];
	}
#endif
	logger(logDebug, "Done")
	return v3;
}

Vec operator +(double a, Vec const & v) {
	logger(logDebug, "a + Vec")
	Vec vr(v.n);
	double *vv = v.values;
	double *vvr = vr.values;
	int n = vr.n;
	for (int i = 0; i < n; i++) {
		vvr[i] = a + vv[i];
	}
	logger(logDebug, "Done")
	return vr;
}

Vec operator -(double a, Vec const & v) {
	logger(logDebug, "a - Vec")
	Vec vr(v.n);
	double *vv = v.values;
	double *vvr = vr.values;
	int n = vr.n;
	for (int i = 0; i < n; i++) {
		vvr[i] = a - vv[i];
	}
	logger(logDebug, "Done")
	return vr;
}

Vec operator *(double a, Vec const & v) {
	logger(logDebug, "a * Vec")
	Vec vr(v.n);
	double *vv = v.values;
	double *vvr = vr.values;
	int n = vr.n;
	for (int i = 0; i < n; i++) {
		vvr[i] = a * vv[i];
	}
	logger(logDebug, "Done")
	return vr;
}

Vec operator /(double a, Vec const & v) {
	logger(logDebug, "a / Vec")
	Vec vr(v.n);
	double *vv = v.values;
	double *vvr = vr.values;
	int n = vr.n;
	for (int i = 0; i < n; i++) {
		vvr[i] = a / vv[i];
	}
	logger(logDebug, "Done")
	return vr;
}

Vec operator +(Vec const & v, double a) {
	logger(logDebug, "Vec + a")
	Vec vr(v.n);
	double *vv = v.values;
	double *vvr = vr.values;
	int n = vr.n;
	for (int i = 0; i < n; i++) {
		vvr[i] = vv[i] + a;
	}
	logger(logDebug, "Done")
	return vr;
}

Vec operator -(Vec const & v, double a) {
	logger(logDebug, "Vec - a")
	Vec vr(v.n);
	double *vv = v.values;
	double *vvr = vr.values;
	int n = vr.n;
	for (int i = 0; i < n; i++) {
		vvr[i] = vv[i] - a;
	}
	logger(logDebug, "Done")
	return vr;
}

Vec operator *(Vec const & v, double a) {
	logger(logDebug, "Vec * a")
	Vec vr(v.n);
	double *vv = v.values;
	double *vvr = vr.values;
	int n = vr.n;
	for (int i = 0; i < n; i++) {
		vvr[i] = vv[i] * a;
	}
	logger(logDebug, "Done")
	return vr;
}

Vec operator /(Vec const & v, double a) {
	logger(logDebug, "Vec / a")
	Vec vr(v.n);
	double *vv = v.values;
	double *vvr = vr.values;
	int n = vr.n;
	for (int i = 0; i < n; i++) {
		vvr[i] = vv[i] / a;
	}
	logger(logDebug, "Done")
	return vr;
}

double operator ^(Vec const & v1, Vec const & v2) {
	double dotProd;
	logger(logDebug, "Vec ^ Vec")
	if (v1.n != v2.n) {
		throw std::invalid_argument(
				"Nonconformant arguments: vectors must have the same size.");
	}
#ifdef WITH_MKL
	dotProd = cblas_ddot(v1.n, v1.values, 1, v2.values, 1);
#else
	double * vv1 = v1.values;
	double * vv2 = v2.values;
	int n = v1.n;
	dotProd = 0;
	for (int i = 0; i < n; i++) {
		dotProd += vv1[i] * vv2[i];
	}
#endif
	logger(logDebug, "Done")
	return dotProd;
}

} /* namespace LinkPred */
