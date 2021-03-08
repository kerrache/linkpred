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

#include "LinkPredConfig.hpp"

#ifdef LINKPRED_WITH_ARMADILLO

#include <linkpred/numerical/numerical.hpp>
#include <linkpred/graphalg/encoders/lem/lem.hpp>
#include <linkpred/utils/log.hpp>
#include <armadillo>
#include <cmath>

namespace LinkPred {

template<typename Network> void LEM<Network>::init() {
	logger(logDebug, "Initialize encoder...")

	logger(logDebug, "Done")
}

template<typename Network> void LEM<Network>::encode() {
	logger(logDebug, "Encoding network...")

	long int nbNodes = net->getNbNodes();
	arma::sp_mat NL(nbNodes, nbNodes);
	for (long int i = 0; i < nbNodes; i++) {
		auto di = net->getDeg(i);
		if (di == 0) {
			throw std::runtime_error("Node with degree 0");
		}
		NL(i, i) = 1;
		for (auto nit = net->neighbBegin(i); nit != net->neighbEnd(i); ++nit) {
			long int j = net->end(*nit);
			auto dj = net->getDeg(j);
			NL(i, j) = -1.0 / (std::sqrt(di * dj));
		}
	}

	arma::vec eigval;
	arma::mat eigvec;

	arma::eigs_sym(eigval, eigvec, NL, dim + 1, "sm", tol);

	codeMap = net->template createNodeMapSP<Vec>();
	for (long int i = 0; i < nbNodes; i++) {
		Vec v(dim);
		for (int k = 0; k < dim; k++) {
			v[k] = eigvec(i, k + 1);
		}
		(*codeMap)[i] = std::move(v);
	}

	logger(logDebug, "Done")
}

#define LEM_CPP
#include "linkpred/instantiations.hpp"
#undef LEM_CPP

} /* namespace LinkPred */

#endif

