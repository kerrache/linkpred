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

#include <linkpred/predictors/ucstpredictor.hpp>
#include "LinkPredConfig.hpp"
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
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UCSTPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::predict(EdgesRandomIteratorT begin,
		EdgesRandomIteratorT end, ScoresRandomIteratorT scores) {

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		*(scores + (it - begin)) = 0;
	}
}

// TODO modify code so that to return exactly l edges when available. For this, either remove randomness or create a random-order non-edge iteratorin the class UNetwork
template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> std::size_t UCSTPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::top(std::size_t l,
		EdgesRandomOutputIteratorT eit, ScoresRandomIteratorT sit) {

	std::size_t k = 0;
#ifdef WITH_MPI
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
#ifdef WITH_MPI
}
#endif
	return k;
}

#define UCSTPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef UCSTPREDICTOR_CPP

}

/* namespace LinkPred */
