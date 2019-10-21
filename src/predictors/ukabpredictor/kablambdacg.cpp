/*
 * This file is part of LinkPred.
 *
 * LinkPred: A high performance library for link prediction in complex simlCalc.getNet()works.
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

#include <linkpred/predictors/ukabpredictor/kablambdacg.hpp>
#include "LinkPredConfig.hpp"
#include "linkpred/utils/randomgen.hpp"
#include <cmath>
#include <iostream>

namespace LinkPred {

void KABLambdaCG::init(double * x, CG::INT n) {
	x[0] = 0.5;
}

double KABLambdaCG::f(double * x, CG::INT n) {

	double lambda = x[0];
	double objp, objn;

	objp = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for reduction (+:objp) if (parallel)
#endif
	for (std::size_t i = 0; i < nbPos; i++) {
		double sum = posSUM[i];
		double nsc = posNSC[i];

		double v = lambda * (sum - nsc) + 1 - sum;
		objp += (v * v) / (2 * nb);
	}

	objn = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for reduction (+:objn) if (parallel)
#endif
	for (std::size_t i = 0; i < nbNeg; i++) {
		double sum = negSUM[i];
		double nsc = negNSC[i];

		double v = lambda * (nsc - sum) + sum;
		objn += (v * v) / (2 * nb);
	}

	return w * objp + (1 - w) * objn;
}

void KABLambdaCG::grad(double * grad_f, double * x, CG::INT n) {
	double lambda = x[0];
	double gp, gn;

	gp = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for reduction (+:gp) if (parallel)
#endif
	for (std::size_t i = 0; i < nbPos; i++) {
		double sum = posSUM[i];
		double nsc = posNSC[i];

		gp += (sum - nsc) * (lambda * (sum - nsc) + 1 - sum) / nb;
	}

	gn = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for reduction (+:gn) if (parallel)
#endif
	for (std::size_t i = 0; i < nbNeg; i++) {
		double sum = negSUM[i];
		double nsc = negNSC[i];

		gn += (nsc - sum) * (lambda * (nsc - sum) + sum) / nb;
	}

	grad_f[0] = w * gp + (1 - w) * gn;
}

double KABLambdaCG::fgrad(double * grad_f, double * x, CG::INT n) {

	grad(grad_f, x, n);
	return f(x, n);
}

void KABLambdaCG::finalize(double * x, CG::INT n) {

	sol = x[0];
	if (sol < 0) {
		sol = 0;
	}
	if (sol > 1) {
		sol = 1;
	}
}

} /* namespace LinkPred */
