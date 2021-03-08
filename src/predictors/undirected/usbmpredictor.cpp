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

#include <linkpred/predictors/undirected/usbmpredictor.hpp>
#include "linkpred/utils/log.hpp"
#include <string>

namespace LinkPred {


template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> struct node_gra * USBMPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::toRGraphNet(
		std::shared_ptr<Network const> net) {

	// Create the header of the graph
	struct node_gra *root = CreateHeaderGraph();
	struct node_gra *prev = root;

	// Adding nodes
	nodesMap.clear();
	nodesMap.reserve(net->getNbNodes());
	for (std::size_t i = 0; i < net->getNbNodes(); i++) {
		prev = CreateNodeGraph(prev, "");
		nodesMap.push_back(prev);
	}

	// Adding edges
	for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it) {
		struct node_gra *n1 = nodesMap[Network::start(*it)];
		struct node_gra *n2 = nodesMap[Network::end(*it)];
		// Add the link
		AddAdjacency(n1, n2, 0, 0, 1, 0);
		AddAdjacency(n2, n1, 0, 0, 1, 0);
	}

	return root;
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void USBMPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::init() {
	logger(logDebug, "Initializing solver...")
	if (rGNet != nullptr) {
		RemoveGraph(rGNet);
	}
	rGNet = toRGraphNet(net);
	logger(logDebug, "Done")
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void USBMPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::learn() {
	logger(logDebug, "Learning...")
	gsl_rng *rng;

	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, seed);

	// Getting scores
	logger(logDebug, "maxIter: " << maxIter)
	scores = LinkScore(rGNet, 0.0, maxIter, rng, 'q');
	gsl_rng_free(rng);
	logger(logDebug, "Done")
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void USBMPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::predict(EdgeRndIt begin,
		EdgeRndIt end, ScoreRndIt oscores) {
	logger(logDebug, "Predicting links...")
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		auto i = Network::start(*it);
		auto j = Network::end(*it);
		*(oscores + (it - begin)) = scores[nodesMap[i]->num][nodesMap[j]->num];
	}
	logger(logDebug, "Done")
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> USBMPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::~USBMPredictor() {
	if (rGNet != nullptr) {
		RemoveGraph(rGNet);
	}

	if (scores != nullptr) {
		free_d_mat(scores, net->getNbNodes());
	}
}

#define USBMPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef USBMPREDICTOR_CPP


} /* namespace LinkPred */
