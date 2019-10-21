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

#include "linkpred/numerical/logistic/logregcg.hpp"
#include "linkpred/utils/randomgen.hpp"

namespace LinkPred {

template class LogRegCG<> ;

template<typename InputRandomIterator, typename OutputRandomIterator> void LogRegCG<
		InputRandomIterator, OutputRandomIterator>::init(double * theta,
		CG::INT n) {

	if (trainingInputsEnd - trainingInputsBegin
			!= trainingOutputsEnd - trainingOutputsBegin) {
		throw std::runtime_error(
				"The size of inputs and outputs are different");
	}
	m = trainingInputsEnd - trainingInputsBegin;
	for (CG::INT i = 0; i < n; i++) {
		theta[i] = rng.getDouble(-1, 1);
	}
}

template<typename InputRandomIterator, typename OutputRandomIterator> double LogRegCG<
		InputRandomIterator, OutputRandomIterator>::f(double * theta,
		CG::INT n) {

	Vec th(n, theta);
	double err = 0;
//#pragma omp parallel for reduction(-:err)
	for (std::size_t i = 0; i < m; i++) {
		double e = 1.0 / (1 + std::exp(-(th ^ *(trainingInputsBegin + i))));
		err -= *(trainingOutputsBegin + i) * std::log(e)
				+ (1 - *(trainingOutputsBegin + i)) * std::log(1 - e);
	}

	double nrm = th.norm();
	double reg = lambda * nrm * nrm;
	return err / m + reg / n;
}

template<typename InputRandomIterator, typename OutputRandomIterator> void LogRegCG<
		InputRandomIterator, OutputRandomIterator>::grad(double * grad_f,
		double * theta, CG::INT n) {
	Vec th(n, theta);
	for (int j = 0; j < n; j++) {
		grad_f[j] = 0;
	}
	for (std::size_t i = 0; i < m; i++) {
		double e = 1.0 / (1 + std::exp(-(th ^ *(trainingInputsBegin + i))));
		for (int j = 0; j < n; j++) {
			grad_f[j] += (e - *(trainingOutputsBegin + i))
					* ((*(trainingInputsBegin + i))[j]);
		}
	}
	for (int j = 0; j < n; j++) {
		grad_f[j] = grad_f[j] / m + 2 * lambda * theta[j] / n;
	}
}

template<typename InputRandomIterator, typename OutputRandomIterator> double LogRegCG<
		InputRandomIterator, OutputRandomIterator>::fgrad(double * grad_f,
		double * theta, CG::INT n) {
	Vec th(n, theta);
	for (int j = 0; j < n; j++) {
		grad_f[j] = 0;
	}
	double err = 0;
	for (std::size_t i = 0; i < m; i++) {
		double e = 1.0 / (1 + std::exp(-(th ^ *(trainingInputsBegin + i))));
		err -= *(trainingOutputsBegin + i) * std::log(e)
				+ (1 - *(trainingOutputsBegin + i)) * std::log(1 - e);
		for (int j = 0; j < n; j++) {
			grad_f[j] += (e - *(trainingOutputsBegin + i))
					* ((*(trainingInputsBegin + i))[j]);
		}
	}
	for (int j = 0; j < n; j++) {
		grad_f[j] = grad_f[j] / m + 2 * lambda * theta[j] / n;
	}
	double nrm = th.norm();
	double reg = lambda * nrm * nrm;
	return err / m + reg / n;
}

}
/* namespace LinkPred */
