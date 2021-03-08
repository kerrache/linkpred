/*
 * This file is part of LinkPred.
 *
 * LinkPred: A high performance library for link prediction in complex networks.
 * Copyright (C) 2021  by Said Kerrache.
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

#include <linkpred/utils/miscutils.hpp>
#include "linkpred/graphalg/encoders/matfact/matfactcg.hpp"
#include "linkpred/utils/randomgen.hpp"
#include <iostream>

#ifdef LINKPRED_WITH_OPENMP
#include <omp.h>
#endif

namespace LinkPred {

MatFactCG::MatFactCG(long int nbNodes,
		std::vector<MatFactPbData> const & pbData, int dim, double lambda,
		long int seed) :
		nbNodes(nbNodes), pbData(pbData), dim(dim), lambda(lambda), seed(seed), rng(
				seed) {

	n = dim * nbNodes; // Number of variables
	m = pbData.size(); // Number of couples
}

void MatFactCG::init(double * x, CG::INT n) {
	double avgTarget = 0;
//#pragma omp parallel for reduction(max:maxD)
	for (auto it = pbData.begin(); it != pbData.end(); ++it) {
		avgTarget += it->target;
	}
	avgTarget /= pbData.size();

	double minCoord = std::sqrt(avgTarget / dim) * 0.5;
	double maxCoord = std::sqrt(avgTarget / dim) * 1.5;

	for (long int i = 0; i < n; i++) {
		x[i] = rng.getDouble(minCoord, maxCoord);
	}
}

double MatFactCG::f(double * x, CG::INT n) {
	// First, compute the MSE
	double mse = 0;
	for (auto it = pbData.begin(); it != pbData.end(); ++it) {
		int offsetI = it->i * dim;
		int offsetJ = it->j * dim;
		double dotProd = 0;
		for (int k = 0; k < dim; k++) {
			dotProd += x[offsetI + k] * x[offsetJ + k];
		}
		double err = dotProd - it->target;
		mse += err * err;
	}
	mse /= m;

	// Second, compute the regularization term
	double reg = 0;
	for (int i = 0; i < n; i++) {
		reg += x[i] * x[i];
	}
	reg /= n;

	double obj = 0.5 * (mse + lambda * reg);
	return obj;
}

void MatFactCG::grad(double * grad_f, double * x, CG::INT n) {
// Initialize to 0
//#pragma omp parallel for
	for (long int i = 0; i < n; i++) {
		grad_f[i] = 0.0;
	}

	// First, compute gradient from the MSE
	for (auto it = pbData.begin(); it != pbData.end(); ++it) {
		int offsetI = it->i * dim;
		int offsetJ = it->j * dim;
		double dotProd = 0;
		for (int k = 0; k < dim; k++) {
			dotProd += x[offsetI + k] * x[offsetJ + k];
		}
		double err = dotProd - it->target;
		for (int k = 0; k < dim; k++) {
			grad_f[offsetI + k] += x[offsetJ + k] * err;
			grad_f[offsetJ + k] += x[offsetI + k] * err;
		}
	}

	for (long int i = 0; i < n; i++) {
		grad_f[i] /= m;
	}

	// Second, compute gradient from the regularization term
	for (int i = 0; i < n; i++) {
		grad_f[i] += lambda * x[i] / n;
	}
}

double MatFactCG::fgrad(double * grad_f, double * x, CG::INT n) {
// Initialize to 0
//#pragma omp parallel for
	for (long int i = 0; i < n; i++) {
		grad_f[i] = 0.0;
	}

	// First, compute gradient from the MSE
	double mse = 0;
	for (auto it = pbData.begin(); it != pbData.end(); ++it) {
		int offsetI = it->i * dim;
		int offsetJ = it->j * dim;
		double dotProd = 0;
		for (int k = 0; k < dim; k++) {
			dotProd += x[offsetI + k] * x[offsetJ + k];
		}
		double err = dotProd - it->target;
		mse += err * err;
		for (int k = 0; k < dim; k++) {
			grad_f[offsetI + k] += x[offsetJ + k] * err;
			grad_f[offsetJ + k] += x[offsetI + k] * err;
		}
	}
	mse /= m;
	for (long int i = 0; i < n; i++) {
		grad_f[i] /= m;
	}

	// Second, compute gradient from the regularization term
	double reg = 0;
	for (int i = 0; i < n; i++) {
		grad_f[i] += lambda * x[i] / n;
		reg += x[i] * x[i];
	}
	reg /= n;

	double obj = 0.5 * (mse + lambda * reg);
	return obj;
}

void MatFactCG::finalize(double * x, CG::INT n) {
}

} /* namespace LinkPred */
