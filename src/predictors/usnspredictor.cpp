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

#include <linkpred/predictors/usnspredictor.hpp>
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
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void USNSPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::predict(EdgesRandomIteratorT begin,
		EdgesRandomIteratorT end, ScoresRandomIteratorT scores) {

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		auto srcNode = NetworkT::start(*it);
		auto endNode = NetworkT::end(*it);
		double ki = net->getDeg(srcNode);
		double kj = net->getDeg(endNode);
		*(scores + (it - begin)) = nsc(srcNode, endNode, ki, kj);
	}
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> double USNSPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::score(EdgeType const & e) {
	auto srcNode = NetworkT::start(e);
	auto endNode = NetworkT::end(e);

	double ki = net->getDeg(srcNode);
	double kj = net->getDeg(endNode);
	return ((1 - lambda) * sum(ki, kj) + lambda * nsc(srcNode, endNode, ki, kj));
}

#define USNSPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef USNSPREDICTOR_CPP

}

/* namespace LinkPred */
