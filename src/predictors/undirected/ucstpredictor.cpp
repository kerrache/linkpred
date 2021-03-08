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

#include <linkpred/predictors/undirected/ucstpredictor.hpp>
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
		typename ScoreRndIt, typename EdgeRndOutIt> void UCSTPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::predict(EdgeRndIt begin,
		EdgeRndIt end, ScoreRndIt scores) {

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		*(scores + (it - begin)) = 0;
	}
}

// TODO modify code so that to return exactly l edges when available. For this, either remove randomness or create a random-order non-edge iteratorin the class UNetwork
template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> std::size_t UCSTPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::top(std::size_t l,
		EdgeRndOutIt eit, ScoreRndIt sit) {

	std::size_t k = 0;
#ifdef LINKPRED_WITH_MPI
	int procID = 0;
	if (distributed) {
		MPI_Comm_rank(comm, &procID);
	}
	if (procID == 0) {
#endif
	double p = (double) l / net->getNbNonEdges();
	for (auto it = net->rndNonEdgesBegin(p, rng.getInt());
			it != net->rndNonEdgesEnd(); ++it) {
		*eit = *it;
		++eit;
		*sit = 0;
		++sit;
		k++;
		if (k >= l) {
			break;
		}
	}
#ifdef LINKPRED_WITH_MPI
}
#endif
	return k;
}

#define UCSTPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef UCSTPREDICTOR_CPP


} /* namespace LinkPred */
