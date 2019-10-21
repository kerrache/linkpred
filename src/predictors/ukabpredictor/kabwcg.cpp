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

#include <linkpred/predictors/ukabpredictor/kabwcg.hpp>
#include "LinkPredConfig.hpp"
#include "linkpred/utils/randomgen.hpp"
#include <cmath>
#include <iostream>

namespace LinkPred {

void KABWCG::init(double * x, CG::INT n) {
	x[0] = 0.5;
}

double KABWCG::f(double * x, CG::INT n) {

	double w = x[0];
	double obj;

	obj = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for reduction (+:obj) if (parallel)
#endif
	for (auto it = posLinksPoNeBegin; it < posLinksPoNeEnd; ++it) {
		double v = 1.0 - (w * it->first + (1 - w) * it->second);
		obj += v * v;
	}
//	std::cout << "#Pos obj: " << obj << std::endl;
//	double posObj = obj;

#ifdef WITH_OPENMP
#pragma omp parallel for reduction (+:obj) if (parallel)
#endif
	for (auto it = negLinksPoNeBegin; it < negLinksPoNeEnd; ++it) {
		double v = (w * it->first + (1 - w) * it->second);
		obj += v * v;
	}

	obj /= (2 * (nbPos + nbNeg));
//	std::cout << "#Neg obj: " << obj - posObj << std::endl;

//	std::cout << "#obj: " << obj << std::endl;
//		std::cout << w << "\t" << obj << std::endl;
//	}
//	exit(0);
	return obj;
}

void KABWCG::grad(double * grad_f, double * x, CG::INT n) {

	double w = x[0];
	double gp = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for reduction (+:gp) if (parallel)
#endif
	for (auto it = posLinksPoNeBegin; it < posLinksPoNeEnd; ++it) {
		gp += (it->second - it->first)
				* (1.0 - (w * it->first + (1 - w) * it->second));
	}
	grad_f[0] = gp;

	double gn = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for reduction (+:gn) if (parallel)
#endif
	for (auto it = negLinksPoNeBegin; it < negLinksPoNeEnd; ++it) {
		gn += (it->first - it->second) * (w * it->first + (1 - w) * it->second);
	}
	grad_f[0] += gn;

	grad_f[0] /= (nbPos + nbNeg);

//	std::cout << "#grad: " << grad_f[0] << std::endl;
}

double KABWCG::fgrad(double * grad_f, double * x, CG::INT n) {

	double w = x[0];
	double obj = 0;
	double gp = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for reduction (+:gp, obj) if (parallel)
#endif
	for (auto it = posLinksPoNeBegin; it < posLinksPoNeEnd; ++it) {
		double v = 1.0 - (w * it->first + (1 - w) * it->second);
		obj += v * v;
		gp += (it->second - it->first) * v;
	}
	grad_f[0] = gp;

	double gn = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for reduction (+:gn, obj) if (parallel)
#endif
	for (auto it = negLinksPoNeBegin; it < negLinksPoNeEnd; ++it) {
		double v = w * it->first + (1 - w) * it->second;
		obj += v * v;
		gn += (it->first - it->second) * (w * it->first + (1 - w) * it->second);
	}
	grad_f[0] += gn;

	obj /= (2 * (nbPos + nbNeg));
	grad_f[0] /= (nbPos + nbNeg);
	return obj;
//	std::cout << "#grad: " << grad_f[0] << std::endl;
}

void KABWCG::finalize(double * x, CG::INT n) {

	sol = x[0];
	if (sol <= 0) {
		sol = 1.0e-8;
	}
	if (sol >= 1) {
		sol = 1 - 1.0e-8;
	}
}

} /* namespace LinkPred */
