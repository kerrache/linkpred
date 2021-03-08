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
#include "linkpred/graphalg/encoders/line/line.hpp"
#include "linkpred/utils/log.hpp"
#include <algorithm>
#include <queue>
namespace LinkPred {

template<typename Network> void LINE<Network>::reconstruct() {
	logger(logDebug, "Reconstruct network...")

	struct Neighbor {
		int nid;
		double weight;
	};

	recNet = std::make_shared<Network>();
	for (long long k = 0; k != nbNodes; k++) {
		recNet->addNode(net->getLabel(k));
	}

	for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it) {
		recNet->addEdge(net->start(*it), net->end(*it));
	}

	std::vector<double> origWDegree(nbNodes); // Weight degree in the original network
	std::fill(origWDegree.begin(), origWDegree.end(), 0);
	if (weightMap) { // Weighted graph
		for (long long i = 0; i != nbNodes; i++) {
			for (auto it = net->neighbBegin(i); it != net->neighbEnd(i);
					++it) {
				origWDegree[i] += weightMap->at(*it);
			}
		}
	} else {  // Un-weighted graph
		for (long long i = 0; i != nbNodes; i++) {
			for (auto it = net->neighbBegin(i); it != net->neighbEnd(i);
					++it) {
				origWDegree[i]++;
			}
		}
	}

	std::map<typename Network::Edge, double> addedEdges;

	for (long long k = 0; k != nbNodes; k++) {
		// Empty the queues
		struct NodeData {
			long long node;
			int depth;
			double weight;
		};
		std::queue < NodeData > q; // BFS queue
		std::map<int, double> vid2weight;

		if (static_cast<int>(net->getDeg(k)) <= threshold) { // We consider only nodes with a small number of edges.

			// This is safe. The default initializer sets value to zero if the key does not exist.
			vid2weight[k] += origWDegree[k] / 5.0; // We divide by 5 instead of 10. In the original code, degree==2*sum_weight

			NodeData nd = { k, 0, origWDegree[k] };
			q.push(nd);

			while (!q.empty()) {
				auto nd = q.front();
				q.pop();

				if (nd.depth != 0) {
					vid2weight[nd.node] += nd.weight;
				}

				if (nd.depth < maxDepth) {
					auto sum = origWDegree[nd.node];
					for (auto it = net->neighbBegin(nd.node);
							it != net->neighbEnd(nd.node); ++it) {
						auto j = net->end(*it);
						double w = 1;
						if (weightMap) {
							w = weightMap->at(*it);
						}
						q.push( { j, nd.depth + 1, nd.weight * w / sum });
					}
				}
			}

			std::vector < Neighbor > rankList;
			rankList.reserve(vid2weight.size());
			for (auto it = vid2weight.begin(); it != vid2weight.end(); it++) {
				if (k != it->first && !net->isEdge(k, it->first)) { // we avoid loops and only add new edges
					rankList.push_back( { it->first, it->second });
				}
			}
			std::sort(rankList.begin(), rankList.end(),
					[](const Neighbor &a, const Neighbor &b) -> bool {
						return a.weight > b.weight;
					});

			long long nbAdded = rankList.size();
			if (nbAdded > threshold) {
				nbAdded = threshold;
			}
			for (int i = 0; i != nbAdded; i++) {
				recNet->addEdge(k, rankList[i].nid);
				addedEdges[net->makeEdge(k, rankList[i].nid)] =
						rankList[i].weight;
			}
		}
	}
	recNet->assemble();

	//  Adding weight for pre-existing edges.
	recWeightMap = recNet->template createEdgeMapSP<double>();
	if (weightMap) {
		for (auto it = weightMap->begin(); it != weightMap->end(); ++it) {
			(*recWeightMap)[it->first] = it->second;
		}
	} else {
		for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it) {
			(*recWeightMap)[*it] = 1;
		}
	}

	// Adding weight for new edges.
	for (auto it = addedEdges.begin(); it != addedEdges.end(); ++it) {
		(*recWeightMap)[it->first] = it->second;
	}

	logger(logDebug, "Done")
}

template<typename Network> void LINE<Network>::init(
		std::shared_ptr<Network const> initNet, WeightMapSP initWeightMap) {
	nbEdges = 2 * initNet->getNbEdges(); // We store both directions of the edge
	edgeFrom.clear();
	edgeFrom.resize(nbEdges);
	edgeTo.clear();
	edgeTo.resize(nbEdges);
	edgeWeight.clear();
	edgeWeight.resize(nbEdges);
	wDegree.clear();
	wDegree.resize(nbNodes);
	std::fill(wDegree.begin(), wDegree.end(), 0);

	if (initWeightMap) { // Weighted graph
		std::size_t k = 0;
		for (NodeID i = 0; i < nbNodes; i++) {
			for (auto it = initNet->neighbBegin(i);
					it != initNet->neighbEnd(i); ++it) {
				auto j = initNet->end(*it);
				edgeFrom[k] = i;
				edgeTo[k] = j;
				auto w = initWeightMap->at(*it);
				edgeWeight[k] = w;
				wDegree[i] += w;
				k++;
			}
		}
	} else {  // Un-weighted graph
		std::size_t k = 0;
		for (NodeID i = 0; i < nbNodes; i++) {
			for (auto it = initNet->neighbBegin(i);
					it != initNet->neighbEnd(i); ++it) {
				auto j = initNet->end(*it);
				edgeFrom[k] = i;
				edgeTo[k] = j;
				wDegree[i]++;
				k++;
			}
			std::fill(edgeWeight.begin(), edgeWeight.end(), 1);
		}
	}

}

template<typename Network> void LINE<Network>::init() {
	logger(logDebug, "Initialize encoder...")

	nbNodes = net->getNbNodes();
	if (enableReconstruct) {
		reconstruct();
		init(recNet, recWeightMap);
	} else {
		init(net, weightMap);
	}

	initNegTable();
	initAliasTable();

	logger(logDebug, "Done")
}

template<typename Network> void LINE<Network>::initAliasTable() {
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

template<typename Network> long long LINE<Network>::sampleEdge() {
	long long k = rng.getUInt(0, nbEdges - 1);
	return rng.getDouble(0, 1) <= prob[k] ? k : alias[k];
}

template<typename Network> void LINE<Network>::initNegTable() {

	negTable.clear();
	negTable.resize(negSize);

	double sum = 0;
	for (NodeID k = 0; k < nbNodes; k++) {
		sum += std::pow(wDegree[k], negSamplingPower);
	}

	double curSum = 0;
	double por = 0;
	int vid = 0;
	for (int k = 0; k < negSize; k++) {
		if ((double) (k + 1) / negSize > por) {
			curSum += std::pow(wDegree[vid], negSamplingPower);
			por = curSum / sum;
			vid++;
		}
		negTable[k] = vid - 1;
	}
}

template<typename Network> void LINE<Network>::update(float *vecU,
		float *vecV, float *vecErr, float lr, const int label) {
	float score = 0; // score = dot(vecU, vecV)
	// TODO enable simd
	for (int d = 0; d < dim; d++) {
		score += vecU[d] * vecV[d];
	}
	float g = (label - sigmoid(score)) * lr;
	for (int d = 0; d < dim; d++) {
		vecErr[d] += g * vecV[d]; // Gradient
	}
	for (int d = 0; d < dim; d++) {
		vecV[d] += score * vecU[d]; // Gradient
	}
}

template<typename Network> void LINE<Network>::encode(int encOrder) {
	vecEmb.clear();
	vecEmb.resize(nbNodes * dim);
	std::generate(vecEmb.begin(), vecEmb.end(), [this]() {
		return (rng.getDouble(0, 1) - 0.5) / dim;
	});
	vecCtx.clear();
	vecCtx.resize(nbNodes * dim);
	std::fill(vecCtx.begin(), vecCtx.end(), 0);

	std::vector<float> vecErr(dim);

	float lr = initLR;
	long long step = 0;
	long long lastStep = 0;

	while (true) {
		if (step > totalSteps)
			break;

		if (step - lastStep > stepInterval) {
			lastStep = step;
			lr = initLR * (1 - step / static_cast<float>(totalSteps + 1));
			// TODO Add decay rate as a data member
			if (lr < initLR * 0.0001)
				lr = initLR * 0.0001;
		}
		step++;

		auto curEdge = sampleEdge();
		auto u = edgeFrom[curEdge];
		auto v = edgeTo[curEdge];

		auto lu = u * dim;
		std::fill(vecErr.begin(), vecErr.end(), 0);

		// Negative sampling
		int target, label;
		for (int k = 0; k < nbNegSamples; k++) {
			if (k == 0) {
				target = v;
				label = 1;
			} else {
				target = negTable[rng.getUInt(0, negSize - 1)];
				label = 0;
			}
			auto lv = target * dim;
			if (encOrder == 1) {
				update(&vecEmb[lu], &vecEmb[lv], &vecErr[0], lr, label);
			} else if (encOrder == 2) {
				update(&vecEmb[lu], &vecCtx[lv], &vecErr[0], lr, label);
			} else {
				throw std::runtime_error(
						"This method only supports order 1 or 2");
			}
		}

		for (int d = 0; d != dim; d++) {
			vecEmb[d + lu] += vecErr[d];
		}
	}

}

template<typename Network> void LINE<Network>::encode() {
	logger(logDebug, "Encoding network...")

	codeMap = net->template createNodeMapSP<Vec>();
	if (order == 1) {
		encode(1);
		for (long long i = 0; i < nbNodes; i++) {
			Vec v(dim);
			for (int d = 0; d < dim; d++) {
				v[d] = vecEmb[i * dim + d];
			}
			(*codeMap)[i] = std::move(v);
		}
	} else if (order == 2) {
		encode(2);
		for (long long i = 0; i < nbNodes; i++) {
			Vec v(dim);
			for (int d = 0; d < dim; d++) {
				v[d] = vecEmb[i * dim + d];
			}
			(*codeMap)[i] = std::move(v);
		}
	} else if (order == 12) {
		int tmpDim = dim; // We split the dimension between the two orders

		dim = (tmpDim + 1) / 2; // Order 1 gets ceil(dim / 2)
		encode(1);
		for (long long i = 0; i < nbNodes; i++) {
			Vec v(dim);
			for (int d = 0; d < dim; d++) {
				v[d] = vecEmb[i * dim + d];
			}
			(*codeMap)[i] = std::move(v);
		}

		dim = tmpDim / 2; // Order 2 gets floor(dim / 2)
		encode(2);
		for (long long i = 0; i < nbNodes; i++) {
			Vec v(dim);
			for (int d = 0; d < dim; d++) {
				v[d] = vecEmb[i * dim + d];
			}
			(*codeMap)[i] = std::move(Vec(codeMap->at(i), v)); // Concatenate first and second order embeddings.
		}

		dim = tmpDim; // Restore dimension
	} else {
		throw std::runtime_error(
				"The parameter order can only take the values 1, 2, and 12");
	}

	logger(logDebug, "Done")
}

#define LINE_CPP
#include "linkpred/instantiations.hpp"
#undef LINE_CPP

} /* namespace LinkPred */
