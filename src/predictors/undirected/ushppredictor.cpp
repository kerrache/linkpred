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

#include <linkpred/predictors/undirected/ushppredictor.hpp>
#include <linkpred/utils/miscutils.hpp>
#include "linkpred/utils/randomgen.hpp"
#include <cmath>
#include <functional>

namespace LinkPred {


template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void USHPPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::init() {
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void USHPPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::createLengthMap() {
	logger(logDebug, "Creating length map...")
	length = net->template createEdgeMapSP<double>();

	for (auto it = net->edgesBegin(); it < net->edgesEnd(); ++it) {
		(*length)[*it] = 1;
	}
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void USHPPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::setLandmarks() {
	logger(logDebug, "Setting landmarks...")
	std::size_t nbNodes = net->getNbNodes();
	std::size_t nbLandmarks = std::max<std::size_t>(1,
			std::round(landmarkRatio * nbNodes));

	logger(logDebug1,
			std::to_string(nbLandmarks) + " landmarks out of "
					+ std::to_string(nbNodes) + " nodes")
	landmarks.clear();
	// Select landmarks according to strategy
	switch (landmarkStrategy) {

	case Random: {
		logger(logDebug1, "Choosing landmarks randomly...")
		if (nbLandmarks < net->getNbNodes()) { // Just to deal with the boundary case
			auto landmarksNodesId = Utils::selectRandom(net->nodesBegin(),
					net->nodesEnd(), nbLandmarks, rng.getInt());

			for (auto it = landmarksNodesId.begin();
					it != landmarksNodesId.end(); ++it) {
				landmarks.insert(it->first);
			}
		} else {

			for (NodeID i = 0; i < net->getNbNodes(); i++) {
				landmarks.insert(i);
			}
		}
		logger(logDebug1, "Done")
		break;
	}

	case Hub: {
		logger(logDebug1, "Choosing hubs as landmarks...")
		std::vector<std::pair<NodeID, std::size_t>> tmpLandmarks;
		Utils::selectTopK<std::pair<NodeID, std::size_t>,
				typename Network::NodeDegIt,
				std::back_insert_iterator<
						std::vector<std::pair<NodeID, std::size_t>>>,
				typename Utils::PairCompRight<NodeID, std::size_t,
						std::greater<std::size_t>> >(net->nodesDegBegin(),
				net->nodesDegEnd(), std::back_inserter(tmpLandmarks),
				nbLandmarks);
		for (auto it = tmpLandmarks.begin(); it != tmpLandmarks.end(); ++it) {
			landmarks.insert(it->first);
		}
		logger(logDebug1, "Done")
		break;
	}

	case IHub: {
		std::vector<std::pair<NodeID, std::size_t>> tmpLandmarks;
		Utils::selectTopK<std::pair<NodeID, std::size_t>,
				typename Network::NodeDegIt,
				std::back_insert_iterator<
						std::vector<std::pair<NodeID, std::size_t>>>,
				typename Utils::PairCompRight<NodeID, std::size_t,
						std::less<std::size_t>> >(net->nodesDegBegin(),
				net->nodesDegEnd(), std::back_inserter(tmpLandmarks),
				nbLandmarks);
		for (auto it = tmpLandmarks.begin(); it != tmpLandmarks.end(); ++it) {
			landmarks.insert(it->first);
		}
		break;
	}

	default:
		throw std::runtime_error("Unknown landmark positioning strategy");
	}

	logger(logDebug, "Done")
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void USHPPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::learn() {
	createLengthMap(); // Create the  edges length map
	if (asp) {
		setLandmarks();
		auto calc = std::make_shared<ASPDistCalculator<Network>>(dijkstra,
				length);
		calc->setLandmarks(landmarks.begin(), landmarks.end());
		distCalc = calc;
	} else {
		distCalc = std::make_shared<ESPDistCalculator<Network>>(dijkstra,
				length, cacheLevel);
	}
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void USHPPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::predict(EdgeRndIt begin,
		EdgeRndIt end, ScoreRndIt scores) {
	logger(logDebug, "Predicting links...")
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		auto srcNode = Network::start(*it);
		auto endNode = Network::end(*it);
		if (net->isEdge(srcNode, endNode)) {
			auto res = distCalc->getIndDist(srcNode, endNode);
			*(scores + (it - begin)) = 1.0 / res.first;
		} else {
			auto res = distCalc->getDist(srcNode, endNode);
			*(scores + (it - begin)) = 1.0 / res.first;
		}
	}
	logger(logDebug, "Done.")
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> double USHPPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::score(Edge const & e) {
	logger(logDebug, "Compute score...")
	auto srcNode = Network::start(e);
	auto endNode = Network::end(e);
	double sc;
	if (net->isEdge(srcNode, endNode)) {
		auto res = distCalc->getIndDist(srcNode, endNode);
		sc = 1.0 / res.first;
	} else {
		auto res = distCalc->getDist(srcNode, endNode);
		sc = 1.0 / res.first;
	}

	logger(logDebug, "Done.")
	return sc;
}

#define USHPPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef USHPPREDICTOR_CPP


} /* namespace LinkPred */
