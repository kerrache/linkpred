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
#include "linkpred/predictors/undirected/ukabpredictor.hpp"
#include <cmath>
#include <stdexcept>
#include <limits>
#ifdef LINKPRED_WITH_OPENMP
#include <omp.h>
#endif
#ifdef LINKPRED_WITH_MPI
#include <mpi.h>
#include <limits>
#include <linkpred/utils/miscutils.hpp>
#endif

namespace LinkPred {


template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UKABPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::init() {

}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UKABPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::createLengthMap() {
	logger(logDebug, "Creating length map...")

	length = net->template createEdgeMapSP<double>();

	for (auto it = net->edgesBegin(); it < net->edgesEnd(); ++it) {
		double w = edgeLength(*it);
		(*length)[*it] = w;
	}
	logger(logDebug, "Done")
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UKABPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::learn() {
	createLengthMap(); // Create the  edges length map
	distCalc = std::make_shared<ESPLDistCalculator<Network>>(dijkstra, length,
			horizLim, cacheLevel);
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UKABPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::predict(EdgeRndIt begin,
		EdgeRndIt end, ScoreRndIt scores) {

	logger(logDebug, "Predicting links...")

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		*(scores + (it - begin)) = score(*it);
	}
	logger(logDebug, "Done")
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> double UKABPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::score(Edge const &e) {

	logger(logDebug, "Predicting links...")

	auto srcNode = Network::start(e);
	auto endNode = Network::end(e);
	auto dist = getDist(srcNode, endNode);
	return score(e, dist);
	logger(logDebug, "Done")
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> double UKABPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::score(Edge const &e, double dist) {

	auto i = Network::start(e);
	auto j = Network::end(e);
	double ki = net->getDeg(i);
	double kj = net->getDeg(j);
	return (sum(ki, kj) + nsc(i, j, ki, kj)) / dist;
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> std::size_t UKABPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::top(std::size_t k,
		EdgeRndOutIt eit, ScoreRndIt sit) {

	std::vector<LMapQueue<Edge, double>> mqs;

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel if (parallel)
	{
#pragma omp critical(initLMapQueueArray)
		{
#endif
	mqs.push_back(LMapQueue<Edge, double>(k));
#ifdef LINKPRED_WITH_OPENMP
		}
	}
#endif

// Computing loop start, end and step in case of distributed processing
#ifdef LINKPRED_WITH_MPI
	int nbProcs;
	int procID;
	if (distributed) {
		MPI_Comm_size(comm, &nbProcs);
		MPI_Comm_rank(comm, &procID);
	} else {
		nbProcs = 1;
		procID = 0;
	}

	auto begin = net->nodesBegin() + (net->getNbNodes() / nbProcs) * procID
			+ std::min(static_cast<std::size_t>(procID),
					net->getNbNodes() % nbProcs);
	auto end = net->nodesBegin() + (net->getNbNodes() / nbProcs) * (procID + 1)
			+ std::min(static_cast<std::size_t>(procID + 1),
					net->getNbNodes() % nbProcs);
#else
	auto begin = net->nodesBegin();
	auto end = net->nodesEnd();
#endif

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel) schedule(dynamic) // We need a dynamic schedule because iteration vary greatly in duration
#endif
	for (auto nit = begin; nit < end; ++nit) {
#ifdef LINKPRED_WITH_OPENMP
		int ind = omp_get_thread_num();
#else
		int ind = 0;
#endif

		auto i = nit->first;
		auto nnzSMap = getFiniteDistMap(i);
		for (auto it = nnzSMap->begin(); it != nnzSMap->end(); ++it) {
			auto j = it->first;
			if (i < j) {
				auto e = net->makeEdge(i, j);
				double sc = score(e, it->second.first);
				mqs[ind].push(e, sc);
			}
		}
	}

#ifdef LINKPRED_WITH_OPENMP
	// Now we merge results obtained by all threads
	if (parallel) {
		LMapQueue<Edge, double>::parMerge(mqs.begin(), mqs.end());
	}
#endif

#ifdef LINKPRED_WITH_MPI
	if (distributed) {
		// Now we merge over all processes
		Utils::merge(mqs[0], comm);
	}
#endif

	for (auto it = mqs[0].begin(); it != mqs[0].end(); ++it) {
		*eit = it->first;
		++eit;
		*sit = it->second;
		++sit;
	}
	return mqs[0].size();
}

#define UKABPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef UKABPREDICTOR_CPP


} /* namespace LinkPred */
