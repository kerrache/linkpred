/*
 * This file is part of LinkPred.
 *
 * LinkPred: A high performance library for link prediction in complex networks.
 * Copyright (C) 2021  by Said Kerrache.
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
#include "linkpred/graphalg/encoders/matfact/matfact.hpp"
#include "linkpred/utils/log.hpp"

namespace LinkPred {

template<typename Network> void MatFact<Network>::init() {
	logger(logDebug, "Initialize encoder...")

	logger(logDebug, "Done")
}

template<typename Network> void MatFact<Network>::encode() {
	logger(logDebug, "Encoding network...")

	std::vector < MatFactPbData > pbData;
	for (auto it = net->rndEdgesBegin(posRatio, rng.getInt());
			it != net->rndEdgesEnd(); ++it) {
		pbData.push_back( { net->start(*it), net->end(*it), 1 }); // Positive example
	}

	for (auto it = net->rndNonEdgesBegin(negRatio, rng.getInt());
			it != net->rndNonEdgesEnd(); ++it) {
		pbData.push_back( { net->start(*it), net->end(*it), 0 }); // Negative example
	}

	long int nbNodes = net->getNbNodes();
	long int n = dim * nbNodes;
	std::vector<double> coords(n);
	CG::cg_stats stats;
	CG::CGDProblem* pb = new MatFactCG(nbNodes, pbData, dim, lambda,
			rng.getInt());
	CG::CGDescent solver(pb);
	solver.cg_descent(coords.data(), n, &stats, nullptr, tol, nullptr, true);

	codeMap = net->template createNodeMapSP<Vec>();
	for (long int i = 0; i < nbNodes; i++) {
		Vec v(dim);
		for (int k = 0; k < dim; k++) {
			v[k] = coords[i * dim + k];
		}
		(*codeMap)[i] = std::move(v);
	}

	logger(logDebug, "Done")
}

#define MATFACT_CPP
#include "linkpred/instantiations.hpp"
#undef MATFACT_CPP

} /* namespace LinkPred */
