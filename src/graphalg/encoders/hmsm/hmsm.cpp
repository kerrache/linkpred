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

#include <linkpred/utils/miscutils.hpp>
#include "linkpred/graphalg/encoders/hmsm/hmsm.hpp"
#include "linkpred/utils/log.hpp"
#include "linkpred/numerical/mds/mds.hpp"

namespace LinkPred {

template<typename Network> void HMSM<Network>::init() {
	logger(logDebug, "Initialize encoder...")
	eucDim = dim - 1;
	logger(logDebug, "Done")
}

template<typename Network> void HMSM<Network>::encode() {
	logger(logDebug, "Encoding network...")

	std::vector<double> coords; // The coordinates of the nodes
	std::vector<double> degProd; // Product of the degrees
	std::vector<double> sqDist; // Squared distances between nodes
	std::vector<double> weight; // Link weight in MDS
	double wPos = 1; // Weight of positive links
	double wNeg = 1; // Weight of negative links
	auto nbNodes = net->getNbNodes();
	auto nbCouples = nbNodes * (nbNodes - 1) / 2;

	double lCoef = std::pow(lProb, -1.0 / alpha) - 1;
	double hCoef = std::pow(hProb, -1.0 / alpha) - 1;
	degProd.resize(nbCouples);
	sqDist.resize(nbCouples);
	weight.resize(nbCouples);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (std::size_t k = 0; k < nbCouples; k++) {
		auto e = net->coupleAtOrd(k);
		double prod = net->getDeg(Network::start(e))
				* net->getDeg(Network::end(e)) + eps;
		degProd[k] = prod;
		if (net->isEdge(e)) {
			sqDist[k] = prod * hCoef * prod * hCoef; // Positive link
			weight[k] = wPos;
		} else {
			sqDist[k] = prod * lCoef * prod * lCoef; // Negative link
			weight[k] = wNeg;
		}
	}

	// Running MDS
	coords.resize(eucDim * nbNodes);
	MDS mds;
	mds.setAlg(MDSAlg::CGMDS);
	mds.setTol(tol);
	mds.setNbRuns(nbMdsRuns);
	mds.solve(nbNodes, nbCouples, sqDist.data(), weight.data(), eucDim,
			coords.data(), rng.getInt(), true);

	Utils::assertNoNaN(coords.cbegin(), coords.cend());

	// center coordinates
	for (int j = 0; j < eucDim; j++) {
		double mean = 0;
		for (NodeID i = 0; i < nbNodes; i++) {
			mean += coords[i * eucDim + j];
		}
		mean /= nbNodes;
		for (NodeID i = 0; i < nbNodes; i++) {
			coords[i * eucDim + j] -= mean;
		}
	}

	codeMap = net->template createNodeMapSP<Vec>();
	for (NodeID i = 0; i < nbNodes; i++) {
		Vec v(dim);
		v[0] = net->getDeg(i);
		for (int j = 0; j < eucDim; j++) {
			v[1 + j] = coords[i * eucDim + j];
		}
		(*codeMap)[i] = std::move(v);
	}

	logger(logDebug, "Done")
}

#define HMSM_CPP
#include "linkpred/instantiations.hpp"
#undef HMSM_CPP

} /* namespace LinkPred */
