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

#include <linkpred/predictors/upatpredictor.hpp>
#include "LinkPredConfig.hpp"
#include "linkpred/core/lmapqueue.hpp"
#ifdef WITH_OPENMP
#include <omp.h>
#endif
#ifdef WITH_MPI
#include <mpi.h>
#include "linkpred/utils/utilities.hpp"
#endif

#include <string>

namespace LinkPred {

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UPATPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::predict(EdgesRandomIteratorT begin,
		EdgesRandomIteratorT end, ScoresRandomIteratorT scores) {
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		auto srcNode = NetworkT::start(*it);
		auto endNode = NetworkT::end(*it);
		*(scores + (it - begin)) = net->getDeg(srcNode) * net->getDeg(endNode);
	}
}

//template<typename NetworkT, typename EdgesRandomIteratorT,
//		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> std::size_t PATPredictor<
//		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
//		EdgesRandomOutputIteratorT>::top(std::size_t l,
//		EdgesRandomOutputIteratorT eit, ScoresRandomIteratorT sit) {
//
//	std::vector<LMapQueue<EdgeType, double>> mqs;
//
//#ifdef WITH_OPENMP
//#pragma omp parallel if (parallel)
//	{
//#pragma omp critical(initLMapQueueArray)
//		{
//#endif
//	mqs.push_back(LMapQueue<EdgeType, double>(l));
//#ifdef WITH_OPENMP
//}
//}
//#endif
//
//	// Computing loop start, end and step in case of distributed processing
//#ifdef WITH_MPI
//	int nbProcs;
//	int procID;
//	MPI_Comm_size(comm, &nbProcs);
//	MPI_Comm_rank(comm, &procID);
//	auto start = net->nonEdgesBegin() + procID;
//	auto end = net->nonEdgesEnd();
//	auto step = nbProcs;
//#else
//	auto start = net->nonEdgesBegin();
//	auto end = net->nonEdgesEnd();
//	auto step = 1;
//#endif
//
//#ifdef WITH_OPENMP
//#pragma omp parallel for if (parallel) schedule(dynamic) // We need a dynamic schedule because iteration vary greatly in duration
//#endif
//	for (auto eit = start; eit < end; eit += step) {
//#ifdef WITH_OPENMP
//		int ind = omp_get_thread_num();
//#else
//		int ind = 0;
//#endif
//		auto i = net->start(*eit);
//		auto j = net->end(*eit);
//		auto sc = net->getDeg(i) * net->getDeg(j);
////		auto mess = "Process " + std::to_string(procID) + " Thread "
////				+ std::to_string(ind) + " : " + std::to_string(i) + "\t"
////				+ std::to_string(j) + "\t" + std::to_string(sc) + "\n";
////		std::cout << mess;
//		mqs[ind].push(*eit, sc);
//	}
//
//#ifdef WITH_OPENMP
//		// Now we merge results obtained by all threads
//		if (parallel) {
//			LMapQueue<EdgeType, double>::parMerge(mqs.begin(), mqs.end());
//		}
//#endif
//
////	for (auto it = mqs[0].begin(); it != mqs[0].end(); ++it) {
////		auto i = net->start(it->first);
////		auto j = net->end(it->first);
////		auto sc = it->second;
////		auto mess = "Process " + std::to_string(procID) + " : "
////				+ std::to_string(i) + "\t" + std::to_string(j) + "\t"
////				+ std::to_string(sc) + "\n";
////		std::cout << mess;
////	}
//#ifdef WITH_MPI
//	// Now we merge over all processes
//if (distributed) {
//	Utilities::merge(mqs[0], comm);
//}
//#endif
//
//	for (auto it = mqs[0].begin(); it != mqs[0].end(); ++it) {
//		*eit = it->first;
//		++eit;
//		*sit = it->second;
//		++sit;
//	}
//	return mqs[0].size();
//}

#define UPATPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef UPATPREDICTOR_CPP

}
/* namespace LinkPred */
