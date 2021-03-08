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
#include "linkpred/core/unetwork/unetwork.hpp"
#include "linkpred/utils/log.hpp"
#include <queue>
#include <omp.h>
#include <string>
#include <sstream>

namespace LinkPred {

template<typename Label, typename NodeID, typename Edge> UNetwork<Label, NodeID,
		Edge>::UNetwork(std::vector<std::pair<Label, Label>> const & edges) {

	logger(logDebug, "Buildling network from edges...")

	for (auto it = edges.begin(); it != edges.end(); ++it) {
		auto istart = addNode(it->first).first;
		auto iend = addNode(it->second).first;
		addEdge(istart, iend);
	}
	assemble();
	logger(logDebug, "Done")
}

template<typename Label, typename NodeID, typename Edge> std::pair<NodeID, bool> UNetwork<
		Label, NodeID, Edge>::addNode(Label const & nodeId) {
	if (assembled) {
		throw std::runtime_error("Cannot add nodes, network already assembled");
	}
	auto res = nodesIdMap.insert(nodeId);
	nbNodes = nodesIdMap.size();
	return res;
}

template<typename Label, typename NodeID, typename Edge> void UNetwork<Label,
		NodeID, Edge>::addEdge(NodeID const & i, NodeID const & j) {
	if (assembled) {
		throw std::runtime_error("Cannot add edges, network already assembled");
	}

	if (i == j) {
		throw std::invalid_argument("Loops are not allowed in the network.");
	}
	edgeSet.insert(makeEdge(i, j));
	edgeSet.insert(makeEdge(j, i));
	nbEdges = edgeSet.size() / 2;
}

template<typename Label, typename NodeID, typename Edge> void UNetwork<Label,
		NodeID, Edge>::assemble() {

	logger(logDebug, "Assembling network...")
	if (assembled) {
		throw std::runtime_error(
				"Cannot assemble network, network already assembled");
	}

	if (nbNodes == 0) {
		throw std::runtime_error(
				"Cannot assemble network, no nodes were added");
	}

	sia.reserve(nbNodes + 1);

	sja.reserve(2 * nbEdges);
	aja.reserve(nbEdges);

	std::size_t prevSI = 0;
	std::size_t prevAI = 0;

	// First edge
	auto eit = edgeSet.begin();
	if (nbEdges > 0) {
		auto e = *eit;
		++eit;
		auto st = start(e);
		auto en = end(e);
		for (NodeID l = 0; l < st; l++) {
			sia.push_back(0);
		}

		sia.push_back(0);
		prevSI = st;

		sja.push_back(e);
		if (st < en) {
			for (NodeID l = 0; l < st; l++) {
				aia.push_back(0);
			}

			aia.push_back(0);
			prevAI = st;

			aja.push_back(e);
		}
	} else {
		sia.push_back(0);
		aia.push_back(0);
	}

	// Remaining edges
	for (std::size_t k = 1; k < 2 * nbEdges; k++) {
		auto e = *eit;
		++eit;
		auto st = start(e);
		auto en = end(e);
		if (st != prevSI) {
			for (NodeID l = prevSI + 1; l < st; l++) {
				sia.push_back(sja.size());
			}

			sia.push_back(sja.size());
			prevSI = st;
		} else if (en == sja[k - 1]) {
			throw std::runtime_error(
					std::string("Cannot assemble network, repeated edge (")
							+ std::to_string(st) + std::string(" , ")
							+ std::to_string(en) + std::string(")"));
		}

		sja.push_back(e);
		if (st < en) {
			if (st != prevAI) {
				for (NodeID l = prevAI + 1; l < st; l++) {
					aia.push_back(aja.size());
				}

				aia.push_back(aja.size());
				prevAI = st;
			}
			aja.push_back(e);
		}
	}
	for (auto i = prevSI + 1; i <= nbNodes; i++) {
		sia.push_back(sja.size());
	}
	for (auto i = prevAI + 1; i <= nbNodes; i++) {
		aia.push_back(aja.size());
	}

	neja.resize(nbEdges + 1);
	if (nbEdges > 0) {
		neja[0] = coupleOrd(aja[0]);
	}
//#pragma omp parallel for
	for (std::size_t i = 1; i < nbEdges; i++) {
		neja[i] = coupleOrd(aja[i]) - coupleOrd(aja[i - 1]) - 1;
	}

	for (std::size_t i = 1; i < nbEdges; i++) {
		neja[i] = neja[i] + neja[i - 1];
	}

	nbNonEdges = nbNodes * (nbNodes - 1) / 2 - nbEdges;
	neja[nbEdges] = nbNonEdges;

	edgeSet.clear();
	assembled = true;

	// Collecting some statistics
	minDeg = std::numeric_limits < std::size_t > ::max();
	maxDeg = 0;
	for (std::size_t i = 0; i < nbNodes; i++) {
		std::size_t deg = sia[i + 1] - sia[i];
		if (minDeg > deg) {
			minDeg = deg;
		}
		if (maxDeg < deg) {
			maxDeg = deg;
		}
	}
	avgDeg = (2.0 * nbEdges) / nbNodes;

	logger(logDebug, "Done")
}

template<typename Label, typename NodeID, typename Edge> void UNetwork<Label,
		NodeID, Edge>::shuffle(long int seed) {

	logger(logDebug, "Shuffling network...")
	if (!assembled) {
		throw std::runtime_error("Network must be assembled before shuffling");
	}

	// Generate random permutation
	auto perm = Utils::getRndPerm(nbNodes, seed);

	// Create new network
	UNetwork<Label, NodeID, Edge> nnet;

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
	this->aia = std::move(nnet.aia);
	this->aja = std::move(nnet.aja);
	this->sia = std::move(nnet.sia);
	this->sja = std::move(nnet.sja);
	this->neja = std::move(nnet.neja);

	logger(logDebug, "Done")
}

template<typename Label, typename NodeID, typename Edge> bool UNetwork<Label,
		NodeID, Edge>::isEdge(Edge const & e) const {

	auto i = start(e);

	if ((i < 0) || (i >= nbNodes)) {
		return false;
	}

	auto l = sia[i];
	auto r = sia[i + 1];
	if (l == r) { // This is added here because the content of sia is unsigned, this way we avoid the problem when l=r=0
		return false;
	}
	r--;

	while (l <= r) {
//		logger(logError, "[" << l << " , " << r << "]")
		auto mid = (l + r) / 2;
		if (e == sja[mid]) {
			return true;
		}
		if (e < sja[mid]) {
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

template<typename Label, typename NodeID, typename Edge> std::shared_ptr<
		UNetwork<Label, NodeID, Edge>> UNetwork<Label, NodeID, Edge>::read(
		std::string fileName, bool ignoreRepetitions, bool ignoreLoops) {
	logger(logDebug, "Reading network from file...")
	std::ifstream in;
	in.open(fileName.c_str(), std::fstream::in);

	if (!in) {
		throw std::runtime_error("Cannot open file: " + fileName);
	}
	auto net = std::make_shared<UNetwork<Label, NodeID, Edge>>();
	Label start, end;
	std::string line;
	while (std::getline(in, line)) {
		Utils::ltrim(line); // Trim from left.
		if (line.length() == 0 || line[0] == '#') {
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

template<typename Label, typename NodeID, typename Edge> std::shared_ptr<
		std::vector<Edge>> UNetwork<Label, NodeID, Edge>::readEdges(
		std::string fileName) const {
	logger(logDebug, "Reading couples from file...")
	std::ifstream in;
	in.open(fileName.c_str(), std::fstream::in);

	if (!in) {
		throw std::runtime_error("Cannot open file: " + fileName);
	}
	auto couples = std::make_shared<std::vector<Edge>>();
	Label start, end;
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

template<typename Label, typename NodeID, typename Edge> std::shared_ptr<
		UNetwork<unsigned int, NodeID, Edge>> UNetwork<Label, NodeID, Edge>::generateERN(
		std::size_t nbNodes, double pr, long int seed) {

	logger(logDebug, "Generating an Erdos-Renyi (random) network...")

	if (nbNodes == 0) {
		throw std::invalid_argument("The number of nodes must not be 0");
	}

	UNetwork<unsigned int, NodeID, Edge> tmpNet;
	auto net = std::make_shared<UNetwork<unsigned int, NodeID, Edge>>();
	for (unsigned int i = 0; i < nbNodes; i++) {
		tmpNet.addNode(i);
		net->addNode(i);
	}
	tmpNet.assemble();

	for (auto it = tmpNet.rndNonEdgesBegin(pr, seed);
			it != tmpNet.rndNonEdgesEnd(); ++it) {
		net->addEdge(start(*it), end(*it));
	}

	net->assemble();
	logger(logDebug, "Done")
	return net;
}

template<typename Label, typename NodeID, typename Edge> std::shared_ptr<
		UNetwork<unsigned int, NodeID, Edge>> UNetwork<Label, NodeID, Edge>::generateRNC(
		std::size_t nbNodes, double pr, long int seed) {

	logger(logDebug, "Generating a connected Erdos-Renyi (random) network...")

	if (nbNodes == 0) {
		throw std::invalid_argument("The number of nodes must not be 0");
	}

	if (pr < 2.0 / nbNodes) {
		throw std::invalid_argument(
				"pr must be >= 2.0 / nbNodes for the graph to be connected");
	}

	RandomGen rng(seed);

	UNetwork<unsigned int, NodeID, Edge> tmpNet;
	auto net = std::make_shared<UNetwork<unsigned int, NodeID, Edge>>();
	for (unsigned int i = 0; i < nbNodes; i++) {
		tmpNet.addNode(i);
		net->addNode(i);
		if (i > 0) {
			unsigned int j = rng.getUInt(0, i - 1);
			tmpNet.addEdge(i, j);
			net->addEdge(i, j);
		}
	}
	tmpNet.assemble();

	for (auto it = tmpNet.rndNonEdgesBegin(pr - 2.0 / nbNodes, rng.getInt());
			it != tmpNet.rndNonEdgesEnd(); ++it) {
		net->addEdge(start(*it), end(*it));
	}

	net->assemble();
	logger(logDebug, "Done")
	return net;
}

template<typename Label, typename NodeID, typename Edge> std::shared_ptr<
		UNetwork<unsigned int, NodeID, Edge>> UNetwork<Label, NodeID, Edge>::generateREG(
		std::size_t h, std::size_t w) {
	logger(logDebug, "Generating an regular grid network...")

	auto net = std::make_shared<UNetwork<unsigned int, NodeID, Edge>>();
	for (unsigned int i = 0; i < h * w; i++) {
		net->addNode(i);
	}

	unsigned int k = 0;
	for (unsigned int i = 0; i < h; i++) {
		for (unsigned int j = 0; j < w; j++) {
			if (i > 0) {
				net->addEdge(k, (i - 1) * w + j);
			}

			if (i < h - 1) {
				net->addEdge(k, (i + 1) * w + j);
			}

			if (j > 0) {
				net->addEdge(k, i * w + j - 1);
			}

			if (j < w - 1) {
				net->addEdge(k, i * w + j + 1);
			}
			k++;
		}
	}

	net->assemble();
	logger(logDebug, "Done")
	return net;
}

template<typename Label, typename NodeID, typename Edge> std::size_t UNetwork<
		Label, NodeID, Edge>::getNbPaths(NodeID const & srcId,
		NodeID const & endId, std::size_t length) const {

	std::queue<std::pair<NodeID, std::size_t>> q;
	q.push(std::make_pair(srcId, 0));
	std::size_t nbPaths = 0;
	while (!q.empty()) {
		auto elem = q.front();
		q.pop();
		if (elem.first == endId && elem.second == length) {
			nbPaths++;
		} else {
			if (elem.second < length) {
				for (auto it = neighbBegin(elem.first);
						it < neighbEnd(elem.first); ++it) {
					q.push(std::make_pair(end(*it), elem.second + 1));
				}
			}
		}
	}
	return nbPaths;
}

template<typename Label, typename NodeID, typename Edge> std::size_t UNetwork<
		Label, NodeID, Edge>::getNbInEdges(std::set<NodeID> const & ns) const {

	std::size_t nb = 0;
	for (auto it1 = ns.cbegin(); it1 != ns.cend(); ++it1) {
		auto it2 = it1;
		for (++it2; it2 != ns.cend(); ++it2) {
			if (isEdge(*it1, *it2)) {
				nb++;
			}
		}
	}
	return nb;
}

template<typename Label, typename NodeID, typename Edge> void UNetwork<Label,
		NodeID, Edge>::write(std::string fileName) const {
	logger(logDebug, "Writing network to file...")
	std::ofstream out(fileName, std::fstream::out);
	if (!out) {
		throw std::runtime_error("Cannot open file: " + fileName);
	}

	for (auto it = aja.begin(); it != aja.end(); ++it) {
		out << getLabel(start(*it)) << "\t" << getLabel(end(*it)) << std::endl;
	}
	out.close();
	logger(logDebug, "Done")
}

template<typename Label, typename NodeID, typename Edge> void UNetwork<Label,
		NodeID, Edge>::print() const {
	logger(logDebug, "Writing network to std::cout...")
	for (auto it = aja.begin(); it != aja.end(); ++it) {
		std::cout << getLabel(start(*it)) << "\t" << getLabel(end(*it))
				<< std::endl;
	}
	logger(logDebug, "Done")
}

#define UNETWORK_CPP
#include "linkpred/instantiations.hpp"
#undef UNETWORK_CPP

} /* namespace LinkPred */
