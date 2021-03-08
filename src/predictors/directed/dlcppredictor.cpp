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
#include "linkpred/predictors/directed/dlcppredictor.hpp"
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
		typename ScoreRndIt, typename EdgeRndOutIt> void DLCPPredictor<
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
		typename ScoreRndIt, typename EdgeRndOutIt> double DLCPPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::score(Edge const & e) {
	auto i = Network::start(e);
	auto j = Network::end(e);
	std::size_t nbcn = 0;
	for (auto it = net->outNeighborsBegin(i); it != net->outNeighborsEnd(i);
			++it) {
		if (net->isEdge(net->end(*it), j)) {
			nbcn++;
		}
	}
	return nbcn + 0.001 * net->getNbPaths(i, j, 3);
}

#define DLCPPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef DLCPPREDICTOR_CPP


} /* namespace LinkPred */
