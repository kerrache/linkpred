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

#include "linkpred/ml/classifiers/ffn/ffn.hpp"
#include "linkpred/utils/miscutils.hpp"
#include "linkpred/utils/log.hpp"
#include <algorithm>
#include <limits>

namespace LinkPred {


template<typename InRndIt, typename OutRndIt,
		typename ScoreRndIt> void FFN<InRndIt,
		OutRndIt, ScoreRndIt>::setAutoArch(int dim) {

	int layerSize = dim;
	Add < mlpack::ann::Linear<> > (layerSize, std::max(layerSize / 2, 2));
	layerSize /= 2;
	while (layerSize > 2) {
		Add<mlpack::ann::SigmoidLayer<> >();
		Add < mlpack::ann::Linear<> > (layerSize, std::max(layerSize / 2, 2));
		layerSize /= 2;
	}
	Add<mlpack::ann::LogSoftMax<> >();
}

template<typename InRndIt, typename OutRndIt,
		typename ScoreRndIt> void FFN<InRndIt,
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

	arma::mat trainLabels(1, m);
	k = 0;
	for (auto it = trOutBegin; it != trOutEnd; ++it) {
		trainLabels(0, k) = (*it) + 1; // As required by mlpack, output must start from 1,
		k++;
	}

	Train(trainData, trainLabels);

	logger(logDebug, "Done")
}

template<typename InRndIt, typename OutRndIt,
		typename ScoreRndIt> void FFN<InRndIt,
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
	arma::mat tmpPred;
	Predict(testData, tmpPred);

	auto it = scoresBegin;
	for (int i = 0; i < k; i++) {
		(*it) = sigmoid(tmpPred(1, i) - tmpPred(0, i)); // Only for binary
		++it;
	}
}

#define FFN_CPP
#include "linkpred/instantiations.hpp"
#undef FFN_CPP

} /* namespace LinkPred */

#endif
