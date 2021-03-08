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

#include "linkpred/ml/classifiers/logistic/logisticregresser.hpp"
#include "linkpred/ml/classifiers/logistic/logregcg.hpp"
#include "linkpred/utils/log.hpp"
#include <limits>

namespace LinkPred {


template<typename InRndIt, typename OutRndIt,
		typename ScoreRndIt> std::vector<Vec> LogisticRegresser<
		InRndIt, OutRndIt, ScoreRndIt>::toOneVec(
		InRndIt begin, InRndIt end) {

	std::vector < Vec > vec;
	if (end > begin) {
		vec.resize(end - begin);
		int k = begin->size();
		int j = 0;
		for (auto it = begin; it < end; ++it) {
			Vec v(k + 1);
			v[0] = 1;
			for (int i = 0; i < k; i++) {
				v[i + 1] = (*it)[i];
			}
			vec[j] = std::move(v);
			j++;
		}
	}
	return vec;
}

template<typename InRndIt, typename OutRndIt,
		typename ScoreRndIt> void LogisticRegresser<
		InRndIt, OutRndIt, ScoreRndIt>::learn(
		InRndIt trInBegin,
		InRndIt trInEnd,
		OutRndIt trOutBegin,
		OutRndIt trOutEnd) {

	logger(logDebug, "Learning...")
	if (trInEnd - trInBegin
			!= trOutEnd - trOutBegin) {
		throw std::runtime_error(
				"The size of inputs and outputs are different");
	}

	m = trInEnd - trInBegin;
	if (m == 0) {
		throw std::runtime_error("No inputs");
	}
	n = trInBegin->size() + 1; // Add one for the constant feature
	theta.resize(n);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = trInBegin; it < trInEnd; ++it) {
		if (n - 1 != static_cast<std::size_t>(it->size())) {
			throw std::runtime_error(
					"Inconsistent number of features in the training set");
		}
	}

	auto input = toOneVec(trInBegin, trInEnd);

	logger(logDebug1, "Scaling all features...")
	maxF.resize(n);
	minF.resize(n);
	minF[0] = 1;
	maxF[0] = 1;

	for (std::size_t j = 1; j < n; j++) {
		maxF[j] = input[0][j];
		minF[j] = input[0][j];
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = input.begin(); it < input.end(); ++it) {
			minF[j] = std::min(minF[j], (*it)[j]);
			maxF[j] = std::max(maxF[j], (*it)[j]);
		}
	}

	double eps = std::numeric_limits<double>::epsilon();
	for (std::size_t j = 1; j < n; j++) {
		double factor = (maxF[j] - minF[j] + eps);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = input.begin(); it < input.end(); ++it) {
			(*it)[j] = ((*it)[j] - minF[j] + eps) / factor;
		}
	}
	logger(logDebug1, "Done")

	CG::cg_stats stats;

	LogRegCG<InRndIt, OutRndIt> pb(input.begin(),
			input.end(), trOutBegin, trOutEnd, lambda,
			seed);

	CG::CGDescent solver(&pb);
	solver.cg_descent(theta.getValues(), n, &stats, nullptr, tol, nullptr,
			true);

//	std::cout << "\tf: " << stats.f << std::endl; /*function value at solution */
//	std::cout << "\tgnorm: " << stats.gnorm << std::endl; /* max abs component of gradient */
//	std::cout << "\titer: " << stats.iter << std::endl; /* number of iterations */
//	std::cout << "\tnfunc: " << stats.nfunc << std::endl; /* number of function evaluations */
//	std::cout << "\tngrad: " << stats.ngrad << std::endl; /* number of gradient evaluations */

	logger(logDebug, "Done")
}

template<typename InRndIt, typename OutRndIt,
		typename ScoreRndIt> void LogisticRegresser<
		InRndIt, OutRndIt, ScoreRndIt>::predict(
		InRndIt testInputsBegin, InRndIt testInputsEnd,
		ScoreRndIt scoresBegin) {

	if (testInputsEnd - testInputsBegin == 0) {
		return;
	}

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = testInputsBegin; it < testInputsEnd; ++it) {
		if (n - 1 != static_cast<std::size_t>(it->size())) {
			throw std::runtime_error(
					"Inconsistent number of features in the training set");
		}
	}

	auto input = toOneVec(testInputsBegin, testInputsEnd);

	double eps = std::numeric_limits<double>::epsilon();
	for (std::size_t j = 1; j < n; j++) {
		double factor = (maxF[j] - minF[j] + eps);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = input.begin(); it < input.end(); ++it) {
			(*it)[j] = ((*it)[j] - minF[j] + eps) / factor;
		}
	}

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = input.begin(); it < input.end(); ++it) {
		*(scoresBegin + (it - input.begin())) = 1.0
				/ (1 + std::exp(-(theta ^ *it)));
	}
}

#define LOGISTICREGRESSER_CPP
#include "linkpred/instantiations.hpp"
#undef LOGISTICREGRESSER_CPP

} /* namespace LinkPred */
