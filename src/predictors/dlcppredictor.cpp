/*
 * This file is part of LinkPred.
 *
 * LinkPred: A high performance library for link prediction in complex networks.
 * Copyright (C) 2017  by LCPd Kerrache.
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
#include "linkpred/predictors/dlcppredictor.hpp"
#include "linkpred/core/lmapqueue.hpp"
#ifdef WITH_OPENMP
#include <omp.h>
#endif
#ifdef WITH_MPI
#include <mpi.h>
#include <limits>
#include "linkpred/utils/utilities.hpp"
#endif

namespace LinkPred {

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void DLCPPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::predict(EdgesRandomIteratorT begin,
		EdgesRandomIteratorT end, ScoresRandomIteratorT scores) {

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
			*(scores + (it - begin)) = score(*it);
		}
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> double DLCPPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::score(EdgeType const & e) {
	auto i = NetworkT::start(e);
	auto j = NetworkT::end(e);
	std::size_t nbcn = 0;
	for (auto it = net->outNeighborsBegin(i); it != net->outNeighborsEnd(i); ++it) {
		if (net->isEdge(net->end(*it), j)) {
			nbcn++;
		}
	}
	return nbcn + 0.001 * net->getNbPaths(i, j, 3);
}

/*template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> std::size_t DLCPPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::top(std::size_t l,
		EdgesRandomOutputIteratorT eit, ScoresRandomIteratorT sit) {

	std::vector<LMapQueue<EdgeType, double>> mqs;

#ifdef WITH_OPENMP
#pragma omp parallel if (parallel)
	{
#pragma omp critical(initLMapQueueArray)
		{
#endif
	mqs.push_back(LMapQueue<EdgeType, double>(l));
#ifdef WITH_OPENMP
}
}
#endif

	// Computing loop start, end and step in case of distributed processing
#ifdef WITH_MPI
	int nbProcs;
	int procID;
	MPI_Comm_size(comm, &nbProcs);
	MPI_Comm_rank(comm, &procID);
	auto start = net->edgesBegin() + procID;
	auto end = net->edgesEnd();
	auto step = nbProcs;
#else
	auto start = net->edgesBegin();
	auto end = net->edgesEnd();
	auto step = 1;
#endif

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel) schedule(dynamic) // We need a dynamic schedule because iteration vary greatly in duration
#endif
	for (auto eit = start; eit < end; eit += step) {
#ifdef WITH_OPENMP
		int ind = omp_get_thread_num();
#else
		int ind = 0;
#endif
		auto i = net->start(*eit);
				auto j = net->end(*eit);
				for (auto nnit = net->neighborsBegin(j); nnit != net->neighborsEnd(j);
						++nnit) {
					auto k = net->end(*nnit);
					if (i < k) {
						auto e = net->makeEdge(i, k);
						if (!net->isEdge(e) && mqs[ind].count(e) == 0) {
							auto sc = score(*it)//net->getNbCommonNeighbors(srcNode, endNode)
							+ epsilon * net->getNbPaths(srcNode, endNode, 3);


									score(*nnit)//net->getNbCommonNeighbors(i, k)
									/ (std::sqrt(
											(double) net->getOutDeg(i) * net->getInDeg(k))
											+ std::numeric_limits<double>::epsilon());
							mqs[ind].push(e, sc);
						}
					}
				}

				j = net->start(*eit);
				i = net->end(*eit);
				for (auto nnit = net->neighborsBegin(j); nnit != net->neighborsEnd(j);
						++nnit) {
					auto k = net->end(*nnit);
					if (i < k) {
						auto e = net->makeEdge(i, k);
						if (!net->isEdge(e) && mqs[ind].count(e) == 0) {
							auto sc = score(*nnit)//net->getNbCommonNeighbors(i, k)
									/ (std::sqrt(
											(double) net->getOutDeg(i) * net->getInDeg(k))
											+ std::numeric_limits<double>::epsilon());
							mqs[ind].push(e, sc);
						}
					}
				}
			}

#ifdef WITH_OPENMP
		// Now we merge results obtained by all threads
		if (parallel) {
			LMapQueue<EdgeType, double>::parMerge(mqs.begin(), mqs.end());
		}
#endif

#ifdef WITH_MPI
	// Now we merge over all processes
	if (distributed) {
		Utilities::merge(mqs[0], comm);
	}
#endif

	for (auto it = mqs[0].begin(); it != mqs[0].end(); ++it) {
		*eit = it->first;
		++eit;
		*sit = it->second;
		++sit;
	}
	return mqs[0].size();
}*/

#define DLCPPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef DLCPPREDICTOR_CPP

}

/* namespace LinkPred */
