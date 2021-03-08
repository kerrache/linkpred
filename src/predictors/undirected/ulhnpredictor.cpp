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

#include <linkpred/predictors/undirected/ulhnpredictor.hpp>
#include "LinkPredConfig.hpp"
#include <limits>
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
		typename ScoreRndIt, typename EdgeRndOutIt> void ULHNPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::predict(EdgeRndIt begin,
		EdgeRndIt end, ScoreRndIt scores) {
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		auto srcNode = Network::start(*it);
		auto endNode = Network::end(*it);
		*(scores + (it - begin)) = (double) net->getNbCommonNeighbors(srcNode,
				endNode)
				/ (net->getDeg(srcNode) * net->getDeg(endNode)
						+ std::numeric_limits<double>::epsilon());
	}
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> std::size_t ULHNPredictor<
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
		for (auto nnit = net->neighbBegin(j); nnit != net->neighbEnd(j);
				++nnit) {
			auto k = net->end(*nnit);
			if (i < k) {
				auto e = net->makeEdge(i, k);
				if (!net->isEdge(e) && mqs[ind].count(e) == 0) {
					auto sc = (double) net->getNbCommonNeighbors(i, k)
							/ (net->getDeg(i) * net->getDeg(k)
									+ std::numeric_limits<double>::epsilon());
					mqs[ind].push(e, sc);
				}
			}
		}

		j = net->start(*eit);
		i = net->end(*eit);
		for (auto nnit = net->neighbBegin(j); nnit != net->neighbEnd(j);
				++nnit) {
			auto k = net->end(*nnit);
			if (i < k) {
				auto e = net->makeEdge(i, k);
				if (!net->isEdge(e) && mqs[ind].count(e) == 0) {
					auto sc = (double) net->getNbCommonNeighbors(i, k)
							/ (net->getDeg(i) * net->getDeg(k)
									+ std::numeric_limits<double>::epsilon());
					mqs[ind].push(e, sc);
				}
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

#define ULHNPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef ULHNPREDICTOR_CPP


} /* namespace LinkPred */
