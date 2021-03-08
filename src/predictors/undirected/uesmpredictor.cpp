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
#include "linkpred/predictors/undirected/uesmpredictor.hpp"
#include "linkpred/utils/log.hpp"

namespace LinkPred {


template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> double UESMPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::score(Edge const &e) {

	auto v1 = encoder->getNodeCode(net->start(e));
	auto v2 = encoder->getNodeCode(net->end(e));

	// Compute similarity
	return simMeasure->sim(v1, v2);
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UESMPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::init() {
	logger(logDebug, "Initializing...")

	// Init encoder
	encoder->init();

	logger(logDebug, "Done")
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UESMPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::learn() {

	// Encode network
	encoder->encode();
}

#define UESMPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef UESMPREDICTOR_CPP


} /* namespace LinkPred */

