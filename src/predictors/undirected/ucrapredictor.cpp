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

#include <linkpred/predictors/undirected/ucrapredictor.hpp>
#include "LinkPredConfig.hpp"
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
		typename ScoreRndIt, typename EdgeRndOutIt> void UCRAPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::predict(EdgeRndIt begin,
		EdgeRndIt end, ScoreRndIt scores) {
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		auto srcNode = Network::start(*it);
		auto endNode = Network::end(*it);

		std::vector<typename Network::Edge> ce; // Common neighbors of both nodes
		net->getCommonNeighbors(srcNode, endNode, std::back_inserter(ce));
		std::vector<typename Network::NodeID> ns;
		for (auto cnit = ce.cbegin(); cnit != ce.cend(); ++cnit) {
			ns.push_back(Network::end(*cnit));
		}
		std::vector < std::size_t > degs; // Internal degree of each of the neighbors
		for (auto it1 = ns.cbegin(); it1 != ns.cend(); ++it1) {
			std::size_t nb = 0;
			for (auto it2 = ns.cbegin(); it2 != ns.cend(); ++it2) {
				if (net->isEdge(*it1, *it2)) {
					nb++;
				}
			}
			degs.push_back(nb);
		}
		double sc = 0;
		for (std::size_t i = 0; i < ns.size(); i++) {
			sc += (double) degs[i] / net->getDeg(ns[i]);
		}
		*(scores + (it - begin)) = sc;
	}
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> std::size_t UCRAPredictor<
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
#pragma omp parallel for if (parallel) schedule(dynamic) // We need a dynamic schedule because iterations vary greatly in duration
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
					std::vector<typename Network::Edge> ce; // Common neighbors of both nodes
					net->getCommonNeighbors(i, k, std::back_inserter(ce));
					std::vector<typename Network::NodeID> ns;
					for (auto it = ce.cbegin(); it != ce.cend(); ++it) {
						ns.push_back(Network::end(*it));
					}
					std::vector < std::size_t > degs; // Internal degree of each of the neighbors
					for (auto it1 = ns.cbegin(); it1 != ns.cend(); ++it1) {
						std::size_t nb = 0;
						for (auto it2 = ns.cbegin(); it2 != ns.cend(); ++it2) {
							if (net->isEdge(*it1, *it2)) {
								nb++;
							}
						}
						degs.push_back(nb);
					}
					double sc = 0;
					for (std::size_t j = 0; j < ns.size(); j++) {
						sc += (double) degs[j] / net->getDeg(ns[j]);
					}
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
					std::vector<typename Network::Edge> ce; // Common neighbors of both nodes
					net->getCommonNeighbors(i, k, std::back_inserter(ce));
					std::vector<typename Network::NodeID> ns;
					for (auto it = ce.cbegin(); it != ce.cend(); ++it) {
						ns.push_back(Network::end(*it));
					}
					std::vector < std::size_t > degs; // Internal degree of each of the neighbors
					for (auto it1 = ns.cbegin(); it1 != ns.cend(); ++it1) {
						std::size_t nb = 0;
						for (auto it2 = ns.cbegin(); it2 != ns.cend(); ++it2) {
							if (net->isEdge(*it1, *it2)) {
								nb++;
							}
						}
						degs.push_back(nb);
					}
					double sc = 0;
					for (std::size_t j = 0; j < ns.size(); j++) {
						sc += (double) degs[j] / net->getDeg(ns[j]);
					}
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

#define UCRAPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef UCRAPREDICTOR_CPP


} /* namespace LinkPred */
