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

#include <linkpred/utils/miscutils.hpp>
#include "linkpred/utils/log.hpp"
#include "linkpred/numerical/plfit/plfit.hpp"
#include "linkpred/utils/randomgen.hpp"
#include "linkpred/numerical/plfit/mt.hpp"
#include <random>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <omp.h>
#include <cmath>
#include <time.h>

namespace LinkPred {
namespace Utils {

std::pair<double, double> plFit(std::vector<std::size_t> const & data) {
	logger(logDebug, "Discrete power law fitting...")
//	print(data, "plfitdata");
	PLFit::plfit_discrete_options_t plfit_discrete_options;
	PLFit::plfit_result_t result;
	PLFit::mt_rng_t rng;
	std::random_device rd;
	PLFit::mt_init(&rng, &rd);
	size_t n = data.size();

	double* ddata = new double[n];
	for (std::size_t i = 0; i < data.size(); i++) {
		ddata[i] = data[i];
	}

	// construct the plfit options
	PLFit::plfit_discrete_options_init(&plfit_discrete_options);
	plfit_discrete_options.rng = &rng;

	// fit the power-law distribution
	plfit_discrete_options.alpha_method = PLFit::PLFIT_LBFGS;
	// Estimate xmin and alpha
	PLFit::plfit_discrete(ddata, n, &plfit_discrete_options,
			&result);

	logger(logDebug1,
			std::to_string(result.alpha) + "\t" + std::to_string(result.xmin)
					+ "\t" + std::to_string(result.L) + "\t"
					+ std::to_string(result.D) + "\t"
					+ std::to_string(result.p))
	delete[] ddata;
	logger(logDebug, "Done")
	return std::make_pair(result.alpha, result.xmin);
}

std::pair<double, double> plFit(std::vector<double> const & data) {
	logger(logDebug, "Continuous power law fitting...")
//	print(data, "plfitdata");
	PLFit::plfit_continuous_options_t plfit_continuous_options;
	PLFit::plfit_result_t result;
	PLFit::mt_rng_t rng;
	std::random_device rd;
	PLFit::mt_init(&rng, &rd);
	size_t n = data.size();

	double* ddata = new double[n];
	for (std::size_t i = 0; i < data.size(); i++) {
		ddata[i] = data[i];
	}

	// construct the plfit options
	plfit_continuous_options_init(&plfit_continuous_options);
	plfit_continuous_options.rng = &rng;

	// fit the power-law distribution
	// Estimate xmin and alpha
	PLFit::plfit_continuous(ddata, n, &plfit_continuous_options,
			&result);
	logger(logDebug1,
			std::to_string(result.alpha) + "\t" + std::to_string(result.xmin)
					+ "\t" + std::to_string(result.L) + "\t"
					+ std::to_string(result.D) + "\t"
					+ std::to_string(result.p))
	delete[] ddata;
	logger(logDebug, "Done")
	return std::make_pair(result.alpha, result.xmin);
}

std::vector<std::size_t> getRndPerm(std::size_t n, long int seed) {
	std::mt19937 g(seed);
	std::vector < std::size_t > perm;
	perm.reserve(n);
	for (std::size_t i = 0; i < n; i++) {
		perm.push_back(i);
	}
	std::shuffle(perm.begin(), perm.end(), g);
	return perm;
}

std::pair<std::size_t, std::size_t> localRange(std::size_t n, int nbProcs,
		int procID) {

	std::size_t begin = (n / nbProcs) * procID
			+ std::min(static_cast<std::size_t>(procID), n % nbProcs);
	std::size_t end = (n / nbProcs) * (procID + 1)
			+ std::min(static_cast<std::size_t>(procID + 1), n % nbProcs);
	return std::make_pair(begin, end);
}

} /* namespace Utils */
} /* namespace LinkPred */

