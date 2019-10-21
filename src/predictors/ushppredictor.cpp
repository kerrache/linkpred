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

#include <linkpred/predictors/ushppredictor.hpp>
#include "linkpred/utils/utilities.hpp"
#include "linkpred/utils/randomgen.hpp"
#include <cmath>
#include <functional>

namespace LinkPred {

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void USHPPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::init() {
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void USHPPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::createLengthMap() {
	logger(logDebug, "Creating length map...")
	length = net->template createEdgeMapSP<double>();

	for (auto it = net->edgesBegin(); it < net->edgesEnd(); ++it) {
		(*length)[*it] = 1;
	}
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void USHPPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::setLandmarks() {
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
			auto landmarksNodesId = Utilities::selectRandom(net->nodesBegin(),
					net->nodesEnd(), nbLandmarks, rng.getInt());

			for (auto it = landmarksNodesId.begin();
					it != landmarksNodesId.end(); ++it) {
				landmarks.insert(it->first);
			}
		} else {

			for (NodeIdType i = 0; i < net->getNbNodes(); i++) {
				landmarks.insert(i);
			}
		}
		logger(logDebug1, "Done")
		break;
	}

	case Hub: {
		logger(logDebug1, "Choosing hubs as landmarks...")
		std::vector<std::pair<NodeIdType, std::size_t>> tmpLandmarks;
		Utilities::selectTopK<std::pair<NodeIdType, std::size_t>,
				typename NetworkT::NodeDegIterator,
				std::back_insert_iterator<
						std::vector<std::pair<NodeIdType, std::size_t>>>,
				typename Utilities::PairCompRight<NodeIdType, std::size_t,
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
		std::vector<std::pair<NodeIdType, std::size_t>> tmpLandmarks;
		Utilities::selectTopK<std::pair<NodeIdType, std::size_t>,
				typename NetworkT::NodeDegIterator,
				std::back_insert_iterator<
						std::vector<std::pair<NodeIdType, std::size_t>>>,
				typename Utilities::PairCompRight<NodeIdType, std::size_t,
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

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void USHPPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::learn() {
	createLengthMap(); // Create the  edges length map
	if (asp) {
		setLandmarks();
		auto calc = std::make_shared < ASPDistCalculator
				< NetworkT >> (dijkstra, length);
		calc->setLandmarks(landmarks.begin(), landmarks.end());
		distCalc = calc;
	} else {
		distCalc = std::make_shared < ESPDistCalculator
				< NetworkT >> (dijkstra, length, cacheLevel);
	}
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void USHPPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::predict(EdgesRandomIteratorT begin,
		EdgesRandomIteratorT end, ScoresRandomIteratorT scores) {
	logger(logDebug, "Predicting links...")
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		auto srcNode = NetworkT::start(*it);
		auto endNode = NetworkT::end(*it);
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

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> double USHPPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::score(EdgeType const & e) {
	logger(logDebug, "Compute score...")
	auto srcNode = NetworkT::start(e);
	auto endNode = NetworkT::end(e);
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

}
/* namespace LinkPred */
