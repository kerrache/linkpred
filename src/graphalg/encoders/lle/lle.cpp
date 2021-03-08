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
#include <linkpred/graphalg/encoders/lle/lle.hpp>
#include <linkpred/utils/log.hpp>
#include <armadillo>

namespace LinkPred {

template<typename Network> void LLE<Network>::init() {
	logger(logDebug, "Initialize encoder...")

	logger(logDebug, "Done")
}

template<typename Network> void LLE<Network>::encode() {
	logger(logDebug, "Encoding network...")

	long int nbNodes = net->getNbNodes();
	arma::sp_mat W(nbNodes, nbNodes);
	for (long int i = 0; i < nbNodes; i++) {
		std::size_t deg = net->getDeg(i);
		W(i, i) = 1;
		if (deg > 0) {
			for (auto nit = net->neighbBegin(i); nit != net->neighbEnd(i);
					++nit) {
				long int j = net->end(*nit);
				W(i, j) = -1.0 / deg;
			}
		}
	}
	arma::sp_mat B = W.t() * W;

	arma::vec eigval;
	arma::mat eigvec;

	arma::eigs_sym(eigval, eigvec, B, dim + 1, "sm", tol);

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

#define LLE_CPP
#include "linkpred/instantiations.hpp"
#undef LLE_CPP

} /* namespace LinkPred */

#endif

