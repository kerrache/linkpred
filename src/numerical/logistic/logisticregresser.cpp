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

#include "LinkPredConfig.hpp"
#include "linkpred/numerical/logistic/logisticregresser.hpp"
#include "linkpred/numerical/logistic/logregcg.hpp"
#include "linkpred/utils/log.hpp"
#include <limits>

namespace LinkPred {

template class LogisticRegresser<> ;

template<typename InputRandomIterator, typename OutputRandomIterator,
		typename ScoreRandomIterator> void LogisticRegresser<
		InputRandomIterator, OutputRandomIterator, ScoreRandomIterator>::learn(
		InputRandomIterator trainingInputsBegin,
		InputRandomIterator trainingInputsEnd,
		OutputRandomIterator trainingOutputsBegin,
		OutputRandomIterator trainingOutputsEnd) {
	logger(logDebug, "Learning...")
	if (trainingInputsEnd - trainingInputsBegin
			!= trainingOutputsEnd - trainingOutputsBegin) {
		throw std::runtime_error(
				"The size of inputs and outputs are different");
	}

	m = trainingInputsEnd - trainingInputsBegin;
	if (m == 0) {
		throw std::runtime_error("No inputs");
	}
	n = trainingInputsBegin->size();
	theta.resize(n);
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = trainingInputsBegin; it < trainingInputsEnd; ++it) {
		if (n != static_cast<std::size_t>(it->size())) {
			throw std::runtime_error(
					"Inconsistent number of features in the training set");
		}
		if ((*it)[0] != 1) {
			throw std::runtime_error("The first feature must be constant 1");
		}
	}

	logger(logDebug1, "Scaling all features...")
	maxF.resize(n);
	minF.resize(n);
	minF[0] = 1;
	maxF[0] = 1;

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (std::size_t j = 1; j < n; j++) {
		maxF[j] = std::numeric_limits<double>::lowest();
		minF[j] = std::numeric_limits<double>::max();
		for (auto iit = trainingInputsBegin; iit < trainingInputsEnd; ++iit) {
			minF[j] = std::min(minF[j], (*iit)[j]);
			maxF[j] = std::max(maxF[j], (*iit)[j]);
		}
	}

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (std::size_t j = 1; j < n; j++) {
		for (auto iit = trainingInputsBegin; iit < trainingInputsEnd; ++iit) {
			(*iit)[j] = ((*iit)[j] - minF[j]
					+ std::numeric_limits<double>::epsilon())
					/ (maxF[j] - minF[j]
							+ std::numeric_limits<double>::epsilon());
		}
	}
	logger(logDebug1, "Done")

	CG::cg_stats stats;

	LogRegCG<InputRandomIterator, OutputRandomIterator> pb(trainingInputsBegin,
			trainingInputsEnd, trainingOutputsBegin, trainingOutputsEnd, lambda,
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

template<typename InputRandomIterator, typename OutputRandomIterator,
		typename ScoreRandomIterator> void LogisticRegresser<
		InputRandomIterator, OutputRandomIterator, ScoreRandomIterator>::predict(
		InputRandomIterator testInputsBegin, InputRandomIterator testInputsEnd,
		ScoreRandomIterator scoresBegin) const {

	if (testInputsEnd - testInputsBegin == 0) {
		return;
	}

	for (auto iit = testInputsBegin; iit < testInputsEnd; ++iit) {
		if (n != static_cast<std::size_t>(iit->size())) {
			throw std::runtime_error(
					"Inconsistent number of features in the training set");
		}
		if ((*iit)[0] != 1) {
			throw std::runtime_error("The first feature must be constant 1");
		}
	}

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (std::size_t j = 1; j < n; j++) {
		for (auto iit = testInputsBegin; iit < testInputsEnd; ++iit) {
			(*iit)[j] = ((*iit)[j] - minF[j]
					+ std::numeric_limits<double>::epsilon())
					/ (maxF[j] - minF[j]
							+ std::numeric_limits<double>::epsilon());
		}
	}

	for (auto iit = testInputsBegin; iit < testInputsEnd; ++iit) {
		*(scoresBegin + (iit - testInputsBegin)) = 1.0
				/ (1 + std::exp(-(theta ^ *iit)));
	}
}

} /* namespace LinkPred */
