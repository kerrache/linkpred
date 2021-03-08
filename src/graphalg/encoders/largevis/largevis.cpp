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

#include "linkpred/graphalg/encoders/largevis/largevis.hpp"
#include "linkpred/utils/log.hpp"
#include <map>
#include <algorithm>

namespace LinkPred {

template<typename Network> void LargeVis<Network>::initAliasTable() {
	alias.clear();
	alias.resize(nbEdges);
	prob.clear();
	prob.resize(nbEdges);

	std::vector<float> normProb(nbEdges);
	std::vector<long long> largeBlock(nbEdges);
	std::vector<long long> smallBlock(nbEdges);

	float sum = 0;
	long long curSmallBlock, curLargeBlock;
	long long numSmallBlock = 0, numLargeBlock = 0;

	for (long long k = 0; k < nbEdges; ++k) {
		sum += edgeWeight[k];

	}

	for (long long k = 0; k < nbEdges; ++k) {
		normProb[k] = edgeWeight[k] * nbEdges / sum;
	}

	for (long long k = nbEdges - 1; k >= 0; --k) {
		if (normProb[k] < 1) {
			smallBlock[numSmallBlock++] = k;
		} else {
			largeBlock[numLargeBlock++] = k;
		}
	}

	while (numSmallBlock && numLargeBlock) {
		curSmallBlock = smallBlock[--numSmallBlock];
		curLargeBlock = largeBlock[--numLargeBlock];
		prob[curSmallBlock] = normProb[curSmallBlock];
		alias[curSmallBlock] = curLargeBlock;
		normProb[curLargeBlock] = normProb[curLargeBlock]
				+ normProb[curSmallBlock] - 1;
		if (normProb[curLargeBlock] < 1) {
			smallBlock[numSmallBlock++] = curLargeBlock;
		} else {
			largeBlock[numLargeBlock++] = curLargeBlock;
		}
	}

	while (numLargeBlock) {
		prob[largeBlock[--numLargeBlock]] = 1;
	}

	while (numSmallBlock) {
		prob[smallBlock[--numSmallBlock]] = 1;
	}
}

template<typename Network> long long LargeVis<Network>::sampleEdge() {
	long long k = rng.getUInt(0, nbEdges - 1);
	return rng.getDouble(0, 1) <= prob[k] ? k : alias[k];
}

template<typename Network> void LargeVis<Network>::initNegTable() {
	long long x, p, i;
	reverse.clear();
	std::vector<long long>(reverse).swap(reverse);
	float sumWeights = 0, dd;
	std::vector<float> weights(nbNodes);

	for (i = 0; i < nbNodes; ++i)
		weights[i] = 0;
	for (x = 0; x < nbNodes; ++x) {
		for (p = head[x]; p >= 0; p = next[p]) {
			weights[x] += edgeWeight[p];
		}
		sumWeights += weights[x] = std::pow(weights[x], negSamplingPower);
	}
	next.clear();
	std::vector<long long>(next).swap(next);
	head.clear();
	negTable.clear();
	negTable.resize(negSize);
	dd = weights[0];
	for (i = x = 0; i < negSize; ++i) {
		negTable[i] = x;
		if (i / (float) negSize > dd / sumWeights && x < nbNodes - 1) {
			dd += weights[++x];
		}
	}
}

template<typename Network> void LargeVis<Network>::encode() {

	logger(logDebug, "Encoding network...")

	vecEmb.clear();
	vecEmb.resize(nbNodes * dim);
	std::generate(vecEmb.begin(), vecEmb.end(),
			[this]() {return (rng.getDouble(0, 1) - 0.5) / dim * 0.0001;});
	long long edgeCountActual = 0;
	long long edgeCount = 0;
	long long lastEdgeCount = 0;
	float lr = initLR;
	float gradClip = 5.0;
	float g, gg;
	std::vector<float> cur(dim);
	std::vector<float> err(dim);

	while (true) {
		if (edgeCount > nbSamples + 2) {
			break;
		}
		if (edgeCount - lastEdgeCount > 10000) {
			edgeCountActual += edgeCount - lastEdgeCount;
			lastEdgeCount = edgeCount;
			lr = initLR * (1 - edgeCountActual / (nbSamples + 1.0));
			if (lr < initLR * 0.0001)
				lr = initLR * 0.0001;
		}
		long long p = sampleEdge();
		long long x = edgeFrom[p];
		long long y = edgeTo[p];
		long long lx = x * dim;
		for (long long i = 0; i < dim; ++i) {
			cur[i] = vecEmb[lx + i], err[i] = 0;
		}
		for (long long i = 0; i < nbNegSamples + 1; ++i) {
			if (i > 0) {
				y = negTable[rng.getUInt(0, negSize - 1)];
				if (y == edgeTo[p])
					continue;
			}
			long long ly = y * dim;
			float f = 0;
			for (long long j = 0; j < dim; ++j) {
				f += (cur[j] - vecEmb[ly + j]) * (cur[j] - vecEmb[ly + j]);
			}
			if (i == 0)
				g = -2 / (1 + f);
			else
				g = 2 * gamma / (1 + f) / (0.1 + f);

			for (long long j = 0; j < dim; ++j) {
				gg = g * (cur[j] - vecEmb[ly + j]);
				if (gg > gradClip)
					gg = gradClip;
				if (gg < -gradClip)
					gg = -gradClip;
				err[j] += gg * lr;

				gg = g * (vecEmb[ly + j] - cur[j]);
				if (gg > gradClip)
					gg = gradClip;
				if (gg < -gradClip)
					gg = -gradClip;
				vecEmb[ly + j] += gg * lr;
			}
		}
		for (long long j = 0; j < dim; ++j) {
			vecEmb[lx + j] += err[j];
		}
		++edgeCount;
	}

	codeMap = net->template createNodeMapSP<Vec>();
	for (long long i = 0; i < nbNodes; i++) {
		Vec v(dim);
		for (int j = 0; j < dim; j++) {
			v[j] = vecEmb[i * dim + j];
		}
		(*codeMap)[i] = std::move(v);
	}

	logger(logDebug, "Done")
}

template<typename Network> void LargeVis<Network>::setNbSamples() {
	if (nbSamples < 0) {
		if (nbNodes < 10000)
			nbSamples = 1000;
		else if (nbNodes < 1000000)
			nbSamples = (nbNodes - 10000) * 9000 / (1000000 - 10000) + 1000;
		else
			nbSamples = nbNodes / 100;
	}
}

template<typename Network> void LargeVis<Network>::init() {

	logger(logDebug, "Initialize encoder...")

	nbNodes = net->getNbNodes();
	nbEdges = 2 * net->getNbEdges(); // We store both directions of the edge
	edgeFrom.resize(nbEdges);
	edgeTo.clear();
	edgeTo.resize(nbEdges);
	edgeWeight.clear();
	edgeWeight.resize(nbEdges);

	if (weightMap) { // Weighted graph
		std::size_t k = 0;
		for (NodeID i = 0; i < nbNodes; i++) {
			for (auto it = net->neighbBegin(i); it != net->neighbEnd(i);
					++it) {
				edgeFrom[k] = net->start(*it);
				edgeTo[k] = net->end(*it);
				edgeWeight[k] = weightMap->at(*it);
				k++;
			}
		}
	} else {  // Un-weighted graph
		std::size_t k = 0;
		for (NodeID i = 0; i < nbNodes; i++) {
			for (auto it = net->neighbBegin(i); it != net->neighbEnd(i);
					++it) {
				edgeFrom[k] = net->start(*it);
				edgeTo[k] = net->end(*it);
				k++;
			}
			std::fill(edgeWeight.begin(), edgeWeight.end(), 1);
		}
	}

	next.clear();
	next.resize(nbEdges);
	std::fill(next.begin(), next.end(), -1);

	head.clear();
	head.resize(nbNodes);
	std::fill(head.begin(), head.end(), -1);

	for (long long p = 0; p < nbEdges; ++p) {
		long long x = edgeFrom[p];
		next[p] = head[x];
		head[x] = p;
	}

	initNegTable();
	initAliasTable();

	logger(logDebug, "Done")
}

#define LARGEVIS_CPP
#include "linkpred/instantiations.hpp"
#undef LARGEVIS_CPP

} /* namespace LinkPred */
