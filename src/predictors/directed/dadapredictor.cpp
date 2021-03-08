/*
 * This file is part of LinkPred.
 *
 * LinkPred: A high performance library for link prediction in complex networks.
 * Copyright (C) 2017  by ADAd Kerrache.
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
#include "linkpred/predictors/directed/dadapredictor.hpp"
#include "linkpred/core/ds/lmapqueue.hpp"
#ifdef LINKPRED_WITH_OPENMP
#include <omp.h>
#endif
#ifdef LINKPRED_WITH_MPI
#include <mpi.h>
#include <limits>
#include "linkpred/utils/miscutils.hpp"
#endif

namespace LinkPred {


template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void DADAPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::predict(EdgeRndIt begin,
		EdgeRndIt end, ScoreRndIt scores) {

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		*(scores + (it - begin)) = score(*it);
	}
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> double DADAPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::score(Edge const & e) {
	auto i = Network::start(e);
	auto j = Network::end(e);
	double sum = 0;
	for (auto it = net->outNeighborsBegin(i); it != net->outNeighborsEnd(i);
			++it) {
		auto k = net->end(*it);
		if (net->isEdge(k, j)) {
			sum += 1.0 / std::log(net->getDeg(*it));
		}
	}
	return sum;
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> std::size_t DADAPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::top(std::size_t l,
		EdgeRndOutIt eit, ScoreRndIt sit) {

	std::vector<LMapQueue<Edge, double>> mqs;

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel if (parallel)
	{
#pragma omp critical(initLMapQueueArray)
		{
#endif
	mqs.push_back(LMapQueue<Edge, double>(l));
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
	auto start = net->edgesBegin() + procID;
	auto end = net->edgesEnd();
	auto step = nbProcs;
#else
	auto start = net->edgesBegin();
	auto end = net->edgesEnd();
	auto step = 1;
#endif

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel) schedule(dynamic) // We need a dynamic schedule because iteration vary greatly in duration
#endif
	for (auto eit = start; eit < end; eit += step) {
#ifdef LINKPRED_WITH_OPENMP
		int ind = omp_get_thread_num();
#else
		int ind = 0;
#endif
		auto i = net->start(*eit);
		auto j = net->end(*eit);
		for (auto nnit = net->inNeighborsBegin(i);
				nnit != net->inNeighborsEnd(i); ++nnit) {
			auto k = net->start(*nnit);
			auto e = net->makeEdge(k, j);
			if (!net->isEdge(e) && mqs[ind].count(e) == 0) {
				auto sc = score(e);
				mqs[ind].push(e, sc);
			}
		}
		for (auto nnit = net->outNeighborsBegin(j);
				nnit != net->outNeighborsEnd(j); ++nnit) {
			auto k = net->end(*nnit);
			auto e = net->makeEdge(i, k);
			if (!net->isEdge(e) && mqs[ind].count(e) == 0) {
				auto sc = score(e);
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
	// Now we merge over all processes
	if (distributed) {
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

#define DADAPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef DADAPREDICTOR_CPP


} /* namespace LinkPred */
