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

#include "linkpred/core/dnetwork.hpp"
#include "linkpred/utils/log.hpp"
#include <queue>
#include <omp.h>
#include <string>
#include <sstream>
#include "linkpred/utils/utilities.hpp"

namespace LinkPred {

template<typename LabelType, typename NodeIdType, typename EdgeType,
		typename Comparator> std::pair<NodeIdType, bool> DNetwork<LabelType,
		NodeIdType, EdgeType, Comparator>::addNode(LabelType const & nodeId) {
	if (assembled) {
		throw std::runtime_error("Cannot add nodes, network already assembled");
	}
	auto res = nodesIdMap.insert(nodeId);
	nbNodes = nodesIdMap.size();
	return res;
}

template<typename LabelType, typename NodeIdType, typename EdgeType,
		typename Comparator> void DNetwork<LabelType, NodeIdType, EdgeType,
		Comparator>::addEdge(NodeIdType const & i, NodeIdType const & j) {
	if (assembled) {
		throw std::runtime_error("Cannot add edges, network already assembled");
	}

	if (i == j) {
		throw std::invalid_argument("Loops are not allowed in the network.");
	}
	outEdgeSet.insert(makeEdge(i, j));
	inEdgeSet.insert(makeEdge(j, i));
	nbEdges = outEdgeSet.size();
}

template<typename LabelType, typename NodeIdType, typename EdgeType,
		typename Comparator> void DNetwork<LabelType, NodeIdType, EdgeType,
		Comparator>::assemble() {

	logger(logDebug, "Assembling network...")
	if (assembled) {
		throw std::runtime_error(
				"Cannot assemble network, network already assembled");
	}

	if (nbNodes == 0) {
		throw std::runtime_error(
				"Cannot assemble network, no nodes were added");
	}

	{
		oia.reserve(nbNodes + 1);
		oja.reserve(nbEdges);

//	std::cout << "--------------------------------------" << std::endl;
//	for (auto it = outEdgeSet.begin(); it != outEdgeSet.end(); ++it) {
//		std::cout << "outEdgeSet: " << start(*it) << "\t" << end(*it)
//				<< std::endl;
//	}
//	std::cout << "--------------------------------------" << std::endl;

		std::size_t prevOI = 0;

		// First edge
		auto eit = outEdgeSet.begin();
		if (nbEdges > 0) {
			auto e = *eit;
			++eit;
			auto st = start(e);
			for (NodeIdType l = 0; l < st; l++) {
				oia.push_back(0);
			}

			oia.push_back(0);
			prevOI = st;
			oja.push_back(e);
		} else {
			oia.push_back(0);
			iia.push_back(0);
		}

		// Remaining edges
		for (std::size_t k = 1; k < nbEdges; k++) {
			auto e = *eit;
			++eit;
			auto st = start(e);
			auto en = end(e);
			if (st != prevOI) {
				for (NodeIdType l = prevOI + 1; l < st; l++) {
					oia.push_back(oja.size());
				}

				oia.push_back(oja.size());
				prevOI = st;
			} else if (en == oja[k - 1]) {
				throw std::runtime_error(
						std::string("Cannot assemble network, repeated edge (")
								+ std::to_string(st) + std::string(" , ")
								+ std::to_string(en) + std::string(")"));
			}

			oja.push_back(e);
		}

		for (auto i = prevOI + 1; i <= nbNodes; i++) {
			oia.push_back(oja.size());
		}

		neja.resize(nbEdges + 1);
		if (nbEdges > 0) {
			neja[0] = coupleOrd(oja[0]);
		}
		for (std::size_t i = 1; i < nbEdges; i++) {
			neja[i] = coupleOrd(oja[i]) - coupleOrd(oja[i - 1]) - 1;
		}

		for (std::size_t i = 1; i < nbEdges; i++) {
			neja[i] = neja[i] + neja[i - 1];
		}

		nbNonEdges = nbNodes * (nbNodes - 1) - nbEdges;
		neja[nbEdges] = nbNonEdges;

		outEdgeSet.clear();
	}

	{
		iia.reserve(nbNodes + 1);
		ija.reserve(nbEdges);

		std::size_t prevII = 0;

		// First edge
		auto eit = inEdgeSet.begin();
		if (nbEdges > 0) {
			auto e = *eit;
			++eit;
			auto st = start(e);
			for (NodeIdType l = 0; l < st; l++) {
				iia.push_back(0);
			}

			iia.push_back(0);
			prevII = st;
			ija.push_back(e);
		} else {
			iia.push_back(0);
			iia.push_back(0);
		}

		// Remaining edges
		for (std::size_t k = 1; k < nbEdges; k++) {
			auto e = *eit;
			++eit;
			auto st = start(e);
			auto en = end(e);
			if (st != prevII) {
				for (NodeIdType l = prevII + 1; l < st; l++) {
					iia.push_back(ija.size());
				}

				iia.push_back(ija.size());
				prevII = st;
			} else if (en == ija[k - 1]) {
				throw std::runtime_error(
						std::string("Cannot assemble network, repeated edge (")
								+ std::to_string(st) + std::string(" , ")
								+ std::to_string(en) + std::string(")"));
			}

			ija.push_back(e);
		}

		for (auto i = prevII + 1; i <= nbNodes; i++) {
			iia.push_back(ija.size());
		}

		inEdgeSet.clear();
	}

	assembled = true;

	// Collecting some statistics
	minDeg = std::numeric_limits<std::size_t>::max();
	maxDeg = 0;
	for (std::size_t i = 0; i < nbNodes; i++) {
		std::size_t deg = oia[i + 1] - oia[i] + iia[i + 1] - iia[i];
		if (minDeg > deg) {
			minDeg = deg;
		}
		if (maxDeg < deg) {
			maxDeg = deg;
		}
	}
	avgDeg = (double) nbEdges / nbNodes;

	// Collecting some statistics
	minOutDeg = std::numeric_limits<std::size_t>::max();
	maxOutDeg = 0;
	for (std::size_t i = 0; i < nbNodes; i++) {
		std::size_t deg = oia[i + 1] - oia[i];
		if (minOutDeg > deg) {
			minOutDeg = deg;
		}
		if (maxOutDeg < deg) {
			maxOutDeg = deg;
		}
	}
	avgOutDeg = (double) nbEdges / nbNodes;

	// Collecting some statistics
	minInDeg = std::numeric_limits<std::size_t>::max();
	maxInDeg = 0;
	for (std::size_t i = 0; i < nbNodes; i++) {
		std::size_t deg = iia[i + 1] - iia[i];
		if (minInDeg > deg) {
			minInDeg = deg;
		}
		if (maxInDeg < deg) {
			maxInDeg = deg;
		}
	}
	avgInDeg = (double) nbEdges / nbNodes;

	logger(logDebug, "Done")
}

template<typename LabelType, typename NodeIdType, typename EdgeType,
		typename Comparator> void DNetwork<LabelType, NodeIdType, EdgeType,
		Comparator>::shuffle(long int seed) {

	logger(logDebug, "Shuffling network...")
	if (!assembled) {
		throw std::runtime_error("Network must be assembled before shuffling");
	}

	// Generate random permutation
	auto perm = Utilities::getRndPerm(nbNodes, seed);

	// Create new network
	DNetwork<LabelType, NodeIdType, EdgeType, Comparator> nnet;

	// Insert nodes
	for (std::size_t i = 0; i < nbNodes; i++) {
		nnet.addNode(nodesIdMap.key(perm[i]));
	}

	// Insert edges
	for (auto it = edgesBegin(); it != edgesEnd(); ++it) {
		nnet.addEdge(nnet.getID(getLabel(start(*it))),
				nnet.getID(getLabel(end(*it))));
	}

	// Assemble network
	nnet.assemble();

	// Copy data members
	this->nodesIdMap = std::move(nnet.nodesIdMap);
	this->iia = std::move(nnet.iia);
	this->ija = std::move(nnet.ija);
	this->oia = std::move(nnet.oia);
	this->oja = std::move(nnet.oja);
	this->neja = std::move(nnet.neja);

	logger(logDebug, "Done")
}

template<typename LabelType, typename NodeIdType, typename EdgeType,
		typename Comparator> bool DNetwork<LabelType, NodeIdType, EdgeType,
		Comparator>::isEdge(EdgeType const & e) const {

	auto i = start(e);

	if ((i < 0) || (i >= nbNodes)) {
		return false;
	}

	auto l = oia[i];
	auto r = oia[i + 1];
	if (l == r) { // This is added here because the content of oia is unsigned, this way we avoid the problem when l=r=0
		return false;
	}
	r--;

	while (l <= r) {
//		logger(logError, "[" << l << " , " << r << "]")
		auto mid = (l + r) / 2;
		if (e == oja[mid]) {
			return true;
		}
		if (e < oja[mid]) {
			if (mid == 0) {
				break;
			}
			r = mid - 1;
		} else {
			l = mid + 1;
		}
	}
	return false;
}

template<typename LabelType, typename NodeIdType, typename EdgeType,
		typename Comparator> std::shared_ptr<
		DNetwork<LabelType, NodeIdType, EdgeType, Comparator>> DNetwork<
		LabelType, NodeIdType, EdgeType, Comparator>::read(std::string fileName,
		bool ignoreRepetitions, bool ignoreLoops) {
	logger(logDebug, "Reading network from file...")
	std::ifstream in;
	in.open(fileName.c_str(), std::fstream::in);

	if (!in) {
		throw std::runtime_error("Cannot open file: " + fileName);
	}
	auto net = std::make_shared<
			DNetwork<LabelType, NodeIdType, EdgeType, Comparator>>();
	LabelType start, end;
	std::string line;
	while (std::getline(in, line)) {
		if (line[0] == '#') {
			continue;
		}
		std::stringstream sline(line);
		sline >> start;
		sline >> end;
		if ((start != end) || !ignoreLoops) {
			auto istart = net->addNode(start).first;
			auto iend = net->addNode(end).first;

			net->addEdge(istart, iend);
		}
	}
	in.close();
	net->assemble();
	logger(logDebug, "Done")
	return net;
}

template<typename LabelType, typename NodeIdType, typename EdgeType,
		typename Comparator> std::shared_ptr<std::vector<EdgeType>> DNetwork<
		LabelType, NodeIdType, EdgeType, Comparator>::readCouples(
		std::string fileName) const {
	logger(logDebug, "Reading couples from file...")
	std::ifstream in;
	in.open(fileName.c_str(), std::fstream::in);

	if (!in) {
		throw std::runtime_error("Cannot open file: " + fileName);
	}
	auto couples = std::make_shared<std::vector<EdgeType>>();
	LabelType start, end;
	while (true) {
		in >> start;
		in >> end;
		if (in.eof()) {
			break;
		}
		auto istart = this->getID(start);
		auto iend = this->getID(end);
		couples->push_back(makeEdge(istart, iend));
	}
	in.close();
	logger(logDebug, "Done")
	return couples;
}

template<typename LabelType, typename NodeIdType, typename EdgeType,
		typename Comparator> void DNetwork<LabelType, NodeIdType, EdgeType,
		Comparator>::write(std::string fileName) const {
	logger(logDebug, "Writing network to file...")
	std::ofstream out(fileName, std::fstream::out);
	if (!out) {
		throw std::runtime_error("Cannot open file: " + fileName);
	}

	for (auto it = oja.begin(); it != oja.end(); ++it) {
		out << getLabel(start(*it)) << "\t" << getLabel(end(*it)) << std::endl;
	}
	out.close();
	logger(logDebug, "Done")
}

template<typename LabelType, typename NodeIdType, typename EdgeType,
		typename Comparator> void DNetwork<LabelType, NodeIdType, EdgeType,
		Comparator>::print() const {
	logger(logDebug, "Writing network to std::cout...")
	for (auto it = oja.begin(); it != oja.end(); ++it) {
		std::cout << getLabel(start(*it)) << "\t" << getLabel(end(*it))
				<< std::endl;
	}
	logger(logDebug, "Done")
}



//by lama
template<typename LabelType, typename NodeIdType, typename EdgeType,
		typename Comparator> std::size_t DNetwork<LabelType, NodeIdType,
		EdgeType, Comparator>::getNbPaths(NodeIdType const & srcId,
		NodeIdType const & endId, std::size_t length) const {

	std::queue<std::pair<NodeIdType, std::size_t>> q;
	q.push(std::make_pair(srcId, 0));
	std::size_t nbPaths = 0;
	while (!q.empty()) {
		auto elem = q.front();
		q.pop();
		if (elem.first == endId && elem.second == length) {
			nbPaths++;
		} else {
			if (elem.second < length) {
				for (auto it = neighborsBegin(elem.first);
						it < neighborsEnd(elem.first); ++it) {
					q.push(std::make_pair(end(*it), elem.second + 1));
				}
			}
		}
	}
	return nbPaths;
}

#define DNETWORK_CPP
#include "linkpred/instantiations.hpp"
#undef DNETWORK_CPP

}
/* namespace LinkPred */
