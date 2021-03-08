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
#include "linkpred/graphalg/encoders/node2vec/node2vec.hpp"
#include "linkpred/utils/log.hpp"
#include <algorithm>

namespace LinkPred {

template<typename Network> void Node2Vec<Network>::update(float *ws,
		float *wt, float *wtCache, float lr, const int label) {
	float score = 0; // score = dot(ws, wt)
	// TODO enable simd
	for (int c = 0; c < dim; c++) {
		score += ws[c] * wt[c];
	}
	score = (label - sigmoid(score)) * lr;
	for (int c = 0; c < dim; c++) {
		wtCache[c] += score * ws[c]; // wt gradient
	}
	for (int c = 0; c < dim; c++) {
		ws[c] += score * wt[c]; // ws gradient
	}
}

template<typename Network> void Node2Vec<Network>::initWalker(int n, int *j,
		float *probs) {
	std::vector<int> smaller, larger;
	for (int i = 0; i < n; i++) {
		if (probs[i] < 1)
			smaller.push_back(i);
		else
			larger.push_back(i);
	}
	while (smaller.size() != 0 && larger.size() != 0) {
		int small = smaller.back();
		smaller.pop_back();
		int large = larger.back();
		larger.pop_back();
		j[small] = large;
		probs[large] += probs[small] - 1;
		if (probs[large] < 1) {
			smaller.push_back(large);
		} else {
			larger.push_back(large);
		}
	}
}

template<typename Network> void Node2Vec<Network>::init() {

	logger(logDebug, "Initialize encoder...")

	nbNodes = net->getNbNodes();
	nbEdges = 2 * net->getNbEdges();
	auto csr = net->getCSR();
	offsets = std::move(csr.first);
	edges = std::move(csr.second);

	totalSteps = nbWalks * nbNodes;
	subSample = 1e-3 * nbNodes * nbWalks * walkLength;

	vecEmb.clear();
	vecEmb.resize(nbNodes * dim);
	std::generate(vecEmb.begin(), vecEmb.end(),
			[this]() {return (rng.getDouble(0, 1) - 0.5) / dim;});
	vecCtx.clear();
	vecCtx.resize(nbNodes * dim);
	std::fill(vecCtx.begin(), vecCtx.end(), 0);

	trainOrder.clear();
	trainOrder.resize(nbNodes);
	for (long long i = 0; i < nbNodes; i++) {
		trainOrder[i] = i;
	}
//	Utils::shuffle(trainOrder.begin(), trainOrder.end(), rng.getInt());
	shuffle(trainOrder, trainOrder.size());

	degrees.clear();
	degrees.resize(nbNodes);
	for (long long i = 0; i < nbNodes; i++) {
		degrees[i] = offsets[i + 1] - offsets[i];
	}

	edgeOffsets.clear();
	edgeOffsets.resize(nbEdges + 1);
	edgeOffsets[0] = 0;
	for (long long i = 0; i < nbEdges; i++) {
		edgeOffsets[i + 1] = edgeOffsets[i] + offsets[edges[i] + 1]
				- offsets[edges[i]];
	}

	n2vQs.clear();
	n2vQs.resize(edgeOffsets[nbEdges]);
	n2vJs.clear();
	n2vJs.resize(edgeOffsets[nbEdges]);
	std::fill(n2vJs.begin(), n2vJs.end(), 0);

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (long long src = 0; src < nbNodes; src++) {
		for (unsigned long long dsti = offsets[src]; dsti < offsets[src + 1];
				dsti++) {
			int dst = edges[dsti];
			double sum = 0;
			int dstDegree = degrees[dst];
			for (unsigned long long dstadji = offsets[dst];
					dstadji < offsets[dst + 1]; dstadji++) {
				long long dstadj = edges[dstadji];
				unsigned long long curidx = edgeOffsets[dsti] + dstadji
						- offsets[dst];
				if (dstadj == src) {
					n2vQs[curidx] = 1 / p;
					sum += 1 / p;
				} else {
					if (hasEdge(dstadj, src)) {
						n2vQs[curidx] = 1;
						sum += 1;
					} else {
						n2vQs[curidx] = 1 / q;
						sum += 1 / q;
					}
				}
			}
//#pragma omp simd
			for (unsigned long long i = edgeOffsets[dsti];
					i < edgeOffsets[dsti] + dstDegree; i++)
				n2vQs[i] *= dstDegree / sum;
			initWalker(dstDegree, &n2vJs[edgeOffsets[dsti]],
					&n2vQs[edgeOffsets[dsti]]);
		}
	}

	// Generating a corpus for negative samples.
	negQs.clear();
	negQs.resize(nbNodes);
	negJs.clear();
	negJs.resize(nbNodes);
	nodeCnts.clear();
	nodeCnts.resize(nbNodes);
	std::fill(negQs.begin(), negQs.end(), 0);

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if(parallel)
#endif
	for (long long i = 0; i < nbNodes * nbWalks; i++) {
		int src = trainOrder[i % nbNodes];
#ifdef LINKPRED_WITH_OPENMP
#pragma omp atomic
#endif
		negQs[src]++;
		if (degrees[src] == 0) {
			continue;
		}
		int lastedgeidx = rng.getSInt(offsets[src], offsets[src + 1] - 1);
		int lastnode = edges[lastedgeidx];
#ifdef LINKPRED_WITH_OPENMP
#pragma omp atomic
#endif
		negQs[lastnode]++;
		for (int j = 2; j < walkLength; j++) {
			if (degrees[lastnode] == 0)
				break;
			lastedgeidx = offsets[lastnode]
					+ walkerDraw(degrees[lastnode],
							&n2vQs[edgeOffsets[lastedgeidx]],
							&n2vJs[edgeOffsets[lastedgeidx]]);
			lastnode = edges[lastedgeidx];
#ifdef LINKPRED_WITH_OPENMP
#pragma omp atomic
#endif
			negQs[lastnode]++;
		}
	}
	for (long long i = 0; i < nbNodes; i++) {
		nodeCnts[i] = negQs[i];
	}
	float sum = 0;
	for (long long i = 0; i < nbNodes; i++) {
		negQs[i] = pow(negQs[i], 0.75f);
		sum += negQs[i];
	}
	for (long long i = 0; i < nbNodes; i++) {
		negQs[i] *= nbNodes / sum;
	}

	initWalker(nbNodes, &negJs[0], &negQs[0]);

	logger(logDebug, "Done")
}

template<typename Network> void Node2Vec<Network>::encode() {
	logger(logDebug, "Encoding network...")
	std::vector<long long> dwrw(walkLength); // we cache one random walk per thread
	std::vector<float> cache(dim); // cache for updating the gradient of a node

	const int startNode = rng.getSInt(0, nbNodes - 1); // This is only useful if there are many threads, otherwise, we can start from node 0
	float lr = initLR;
	long long step = 0;
	long long lastStep = 0;

	while (true) {
		if (step > totalSteps)
			break;
		if (step - lastStep > stepInterval) {
			lr = initLR * (1 - step / static_cast<float>(totalSteps + 1)); // linear LR decay
			if (lr < initLR * 0.0001) {
				lr = initLR * 0.0001;
			}
		}
		step++;

		dwrw[0] = trainOrder[(step + startNode) % nbNodes];
		if (degrees[dwrw[0]] == 0) {
			continue;
		}
		unsigned long long lastedgeidx = offsets[dwrw[0]]
				+ rng.getSInt(0, offsets[dwrw[0] + 1] - offsets[dwrw[0]] - 1);
		dwrw[1] = edges[lastedgeidx];
		for (int i = 2; i < walkLength; i++) {
			int lastNode = dwrw[i - 1];
			if (degrees[lastNode] == 0) {
				dwrw[i] = -2;
				break;
			}
			lastedgeidx = offsets[lastNode]
					+ walkerDraw(degrees[lastNode],
							&n2vQs[edgeOffsets[lastedgeidx]],
							&n2vJs[edgeOffsets[lastedgeidx]]);
			dwrw[i] = edges[lastedgeidx];
		}

		for (int dwi = 0; dwi < walkLength; dwi++) {
			int b = rng.getSInt(0, windowSize - 1); // Sub-sample window size
			if (dwrw[dwi] < 0) {
				break;
			}
			size_t n1 = dwrw[dwi];
			if ((std::sqrt(nodeCnts[n1] / subSample) + 1) * subSample
					/ nodeCnts[n1] < rng.getDouble(0, 1)) { // Randomly sub-sample frequent nodes
				continue;
			}

			for (int dwj = std::max(0, dwi - windowSize + b);
					dwj < std::min(dwi + windowSize - b + 1, walkLength);
					dwj++) {
				if (dwi == dwj) {
					continue;
				}
				if (dwrw[dwj] < 0) {
					break;
				}
				size_t n2 = dwrw[dwj];

				std::fill(cache.begin(), cache.end(), 0); // clear cache
				update(&vecCtx[n1 * dim], &vecEmb[n2 * dim], &cache[0], lr, 1);
				for (int i = 0; i < nbNegSamples; i++) {
					size_t neg = walkerDraw(nbNodes, &negQs[0], &negJs[0]);
					while (neg == n2) {
						neg = walkerDraw(nbNodes, &negQs[0], &negJs[0]);
					}
					update(&vecCtx[neg * dim], &vecEmb[n2 * dim], &cache[0], lr,
							0);
				}

				for (int c = 0; c < dim; c++) {
					vecEmb[n2 * dim + c] += cache[c];
				}
			}
		}
	}

	codeMap = net->template createNodeMapSP<Vec>();
	for (long long i = 0; i < nbNodes; i++) {
		Vec v(dim);
		for (int d = 0; d < dim; d++) {
			v[d] = vecEmb[i * dim + d];
		}
		(*codeMap)[i] = std::move(v);
	}
	logger(logDebug, "Done")
}

#define NODE2VEC_CPP
#include "linkpred/instantiations.hpp"
#undef NODE2VEC_CPP

} /* namespace LinkPred */
