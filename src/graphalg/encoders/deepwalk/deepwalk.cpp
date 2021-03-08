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
#include "linkpred/graphalg/encoders/deepwalk/deepwalk.hpp"
#include "linkpred/utils/log.hpp"
#include <algorithm>
#include <numeric>

namespace LinkPred {

template<typename Network> void DeepWalk<Network>::update(float *ws,
		float *wt, float *wtCache, float lr, const int label) {
	float score = 0; // score = dot(ws, wt)
	// TODO Enable simd
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

template<typename Network> void DeepWalk<Network>::estimatePrRW(
		std::vector<float> &outputs, int nbSamples, float alpha) {

	for (long long k = 0; k < nbNodes; k++) {
		outputs[k] = 1;
	}

	for (int i = 0; i < nbSamples; i++) {
		int currentNode = rng.getUInt(0, nbNodes - 1);
		outputs[currentNode]++;
		while (rng.getDouble(0, 1) < alpha) { // kill with probability 1-alpha
			currentNode = sampleNeighbor(currentNode);
			if (currentNode == -1)
				break;
			outputs[currentNode]++;
		}
	}
}

template<typename Network> void DeepWalk<Network>::initHsm(
		std::vector<float> &probs) {

	std::vector < std::size_t > idx(nbNodes); // index array for vertices
	std::iota(idx.begin(), idx.end(), 0);
	sort(idx.begin(), idx.end(), // we need to sort the index array as in probs
			[&probs](size_t i1, size_t i2) {
				return probs[i1] > probs[i2];
			});
	std::vector<float> count(nbNodes * 2 + 1);
	std::vector < uint8_t > binary(nbNodes * 2 + 1);
	std::vector<int> parentNode(nbNodes * 2 + 1);
	for (int k = 0; k < nbNodes; k++) {
		count[k] = probs[idx[k]];
	}
	for (int k = nbNodes; k < nbNodes * 2; k++) {
		count[k] = 1e25;
	}
	int pos1 = nbNodes - 1;
	int pos2 = nbNodes;
	int min1i, min2i; // relentlessly copied from word2vec
	for (int k = 0; k < nbNodes - 1; k++) {
		// First, find two smallest nodes 'min1, min2'
		if (pos1 >= 0) {
			if (count[pos1] < count[pos2]) {
				min1i = pos1;
				pos1--;
			} else {
				min1i = pos2;
				pos2++;
			}
		} else {
			min1i = pos2;
			pos2++;
		}
		if (pos1 >= 0) {
			if (count[pos1] < count[pos2]) {
				min2i = pos1;
				pos1--;
			} else {
				min2i = pos2;
				pos2++;
			}
		} else {
			min2i = pos2;
			pos2++;
		}
		count[nbNodes + k] = count[min1i] + count[min2i];
		parentNode[min1i] = nbNodes + k;
		parentNode[min2i] = nbNodes + k;
		binary[min2i] = 1;
	}

	hsmIndPtrs.clear();
	hsmIndPtrs.resize(nbNodes + 1);
	size_t totalLen = 0;
	for (int k = 0; k < nbNodes; k++) {
		int b = k;
		int i = 0;
		while (true) {
			totalLen++;
			i++;
			b = parentNode[b];
			if (b == nbNodes * 2 - 2)
				break;
		}
		hsmIndPtrs[idx[k]] = i;
	}
	hsmIndPtrs[nbNodes] = totalLen;

	for (int k = nbNodes - 1; k >= 0; k--) {
		hsmIndPtrs[k] = hsmIndPtrs[k + 1] - hsmIndPtrs[k];
	}

	hsmPtrs.clear();
	hsmPtrs.resize(totalLen + 1);
	hsmCodes.clear();
	hsmCodes.resize(totalLen + 1);
	std::fill(hsmCodes.begin(), hsmCodes.end(), 0);
	std::vector<int> point(MaxCodeLength);
	std::vector < uint8_t > code(MaxCodeLength);

	for (long long k = 0; k < nbNodes; k++) {
		long long b = k;
		long long i = 0;
		while (true) {
			code[i] = binary[b];
			point[i] = b;
			i++;
			b = parentNode[b];
			if (b == nbNodes * 2 - 2)
				break;
		}
		int ida = idx[k];
		size_t curptr = hsmIndPtrs[ida];
		for (b = 0; b < i; b++) {
			hsmCodes[curptr + i - b - 1] = code[b];
			hsmPtrs[curptr + i - b] = point[b] - nbNodes;
		}
	}

	for (int k = 0; k < nbNodes; k++) {
		hsmPtrs[hsmIndPtrs[idx[k]]] = nbNodes - 2;
	}
}

template<typename Network> void DeepWalk<Network>::init() {

	logger(logDebug, "Initialize encoder...")

	nbNodes = net->getNbNodes();
	nbEdges = 2 * net->getNbEdges();
	auto csr = net->getCSR();
	offsets = std::move(csr.first);
	edges = std::move(csr.second);

	totalSteps = nbWalks * nbNodes;

	vecEmb.clear();
	vecEmb.resize(nbNodes * dim);
	std::generate(vecEmb.begin(), vecEmb.end(), [this]() {
		return (rng.getDouble(0, 1) - 0.5) / dim;
	});
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

	hsmWeights.clear();
	hsmWeights.resize(nbNodes);
	estimatePrRW(hsmWeights, nbNodes * nbNodeWalks, alpha);
	initHsm(hsmWeights);

	logger(logDebug, "Done")
}

template<typename Network> long long DeepWalk<Network>::sampleNeighbor(
		long long i) {
	if (degrees[i] == 0) {
		return -1;
	}
	auto j = rng.getUInt(0, degrees[i] - 1);
	return net->end(*(net->neighbBegin(i) + j));
}

template<typename Network> void DeepWalk<Network>::encode() {
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
		for (int i = 1; i < walkLength; i++) {
			dwrw[i] = sampleNeighbor(dwrw[i - 1]);
			if (dwrw[i] == -1) {
				break;
			}
		}

		for (int dwi = 0; dwi < walkLength; dwi++) {
			int b = rng.getSInt(0, windowSize - 1); // Sub-sample window size
			size_t n1 = dwrw[dwi];
			if (n1 < 0) {
				break;
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

				for (size_t hsi = hsmIndPtrs[n1]; hsi < hsmIndPtrs[n1 + 1];
						hsi++) {
					size_t tou = hsmPtrs[hsi]; // pointer at level hsi
					int label = hsmCodes[hsi]; // label at level hsi - hsm_indptrs[n1]
					update(&vecCtx[tou * dim], &vecEmb[n2 * dim], &cache[0], lr,
							label);
				}

				for (int d = 0; d < dim; d++) {
					vecEmb[n2 * dim + d] += cache[d];
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

#define DEEPWALK_CPP
#include "linkpred/instantiations.hpp"
#undef DEEPWALK_CPP

} /* namespace LinkPred */
