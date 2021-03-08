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

#include <linkpred/predictors/undirected/ueclpredictor.hpp>
#include "LinkPredConfig.hpp"
#include "linkpred/utils/log.hpp"

namespace LinkPred {


template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> double UECLPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::score(Edge const &e) {
	std::vector<Vec> input;
	std::vector<double> output;
	input.emplace_back(encoder->getEdgeCode(e));
	output.resize(1); // Allocate memory

	// Classify
	classifier->predict(input.begin(), input.end(), output.begin());
	return output[0];
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UECLPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::init() {
	logger(logDebug, "Initializing...")

	// Init encoder
	encoder->init();

	logger(logDebug, "Done")
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UECLPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::learn() {

	// Encode network
	encoder->encode();

	// Prepare training data
	std::vector<Vec> input;
	std::vector<bool> output;

	// First positive edges
	for (auto it = net->rndEdgesBegin(posRatio, rng.getInt());
			it != net->rndEdgesEnd(); ++it) {
		input.push_back(encoder->getEdgeCode(*it));
		output.push_back(1); // Positive example
	}

	// Now negative edges
	for (auto it = net->rndNonEdgesBegin(negRatio, rng.getInt());
			it != net->rndNonEdgesEnd(); ++it) {
		input.push_back(encoder->getEdgeCode(*it));
		output.push_back(0); // Negative example
	}

	// Train the classifier
	classifier->learn(input.begin(), input.end(), output.begin(), output.end());
}

#define UECLPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef UECLPREDICTOR_CPP


} /* namespace LinkPred */

