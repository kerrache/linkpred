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

#ifdef LINKPRED_WITH_MLPACK

#include "linkpred/ml/classifiers/naivebayes/naivebayes.hpp"
#include "linkpred/utils/log.hpp"
#include <algorithm>
#include <limits>

namespace LinkPred {


template<typename InRndIt, typename OutRndIt,
		typename ScoreRndIt> void NaiveBayes<InRndIt,
		OutRndIt, ScoreRndIt>::learn(
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
	n = trInBegin->size();
	for (auto it = trInBegin; it < trInEnd; ++it) {
		if (n != static_cast<std::size_t>(it->size())) {
			throw std::runtime_error(
					"Inconsistent number of features in the training set");
		}
	}

	arma::mat trainData(n, m);
	int k = 0;
	for (auto it = trInBegin; it != trInEnd; ++it) {
		for (std::size_t i = 0; i < n; i++) {
			trainData(i, k) = (*it)[i];
		}
		k++;
	}

	arma::Row<unsigned long int> trainLabels(m);
	k = 0;
	for (auto it = trOutBegin; it != trOutEnd; ++it) {
		trainLabels[k] = (*it); // As required by mlpack, output must start from 0,
		k++;
	}

	Train(trainData, trainLabels, 2);

	logger(logDebug, "Done")
}

template<typename InRndIt, typename OutRndIt,
		typename ScoreRndIt> void NaiveBayes<InRndIt,
		OutRndIt, ScoreRndIt>::predict(
		InRndIt testInputsBegin, InRndIt testInputsEnd,
		ScoreRndIt scoresBegin) {

	if (testInputsEnd - testInputsBegin == 0) {
		return;
	}

	for (auto it = testInputsBegin; it < testInputsEnd; ++it) {
		if (n != static_cast<std::size_t>(it->size())) {
			throw std::runtime_error(
					"Inconsistent number of features in the training set");
		}
	}

	arma::mat testData(n, testInputsEnd - testInputsBegin);
	int k = 0;
	for (auto it = testInputsBegin; it != testInputsEnd; ++it) {
		for (int i = 0; i < it->size(); i++) {
			testData(i, k) = (*it)[i];
		}
		k++;
	}

	arma::Row<unsigned long int> pred;
	arma::mat prob;

	Classify(testData, pred, prob);

	auto it = scoresBegin;
	for (int i = 0; i < k; i++) {
		(*it) = prob(1, i); // We take the probability of the positive class.
		++it;
	}
}

#define NAIVEBAYES_CPP
#include "linkpred/instantiations.hpp"
#undef NAIVEBAYES_CPP

} /* namespace LinkPred */

#endif
