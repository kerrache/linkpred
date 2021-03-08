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

#include <linkpred.hpp>
#include <set>
#include <iostream>
using namespace LinkPred;

/**
 * @brief A class that collects nodes during traversal up to a certain limit.
 * @tparam Network The network type.
 */
template<typename Network = DNetwork<>> class CollectorLim {
protected:
	std::queue<typename Network::NodeID> visited; /**< To store visited nodes. */
	std::size_t limit = 1; /**< Maximum number of nodes to be collected. */

public:

	/**
	 * Constructor.
	 */
	CollectorLim(std::size_t limit) :
			limit(limit) {
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	CollectorLim(CollectorLim const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	CollectorLim & operator =(CollectorLim const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	CollectorLim(CollectorLim && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	CollectorLim & operator =(CollectorLim && that) = default;

	/**
	 * Node processing.
	 * @param i The node's ID.
	 */
	bool process(typename Network::NodeID const & i) {
		visited.push(i);
		return visited.size() < limit;
	}

	/**
	 * @return The visited nodes.
	 */
	const std::queue<typename Network::NodeID>& getVisited() const {
		return visited;
	}

	/**
	 * Destructor.
	 */
	virtual ~CollectorLim() = default;

};

int main(int argc, char*argv[]) {
	if (argc != 4) {
		std::cerr << "Bad arguments\nUsage: " << argv[0]
				<< " netFileName srcNode limit\n";
		exit(1);
	}
	auto net = DNetwork<>::read(std::string(argv[1]), false, true);
	std::string srcNode = argv[2];
	std::size_t limit = std::atol(argv[3]);

	// BFS
	BFS<DNetwork<>, CollectorLim<>> bfs(net);
	CollectorLim<> col(limit);
	bfs.traverse(net->getID(srcNode), col);
	auto visited = col.getVisited();
	std::set<typename DNetwork<>::NodeID> visitedSet;
	while (!visited.empty()) {
		auto i = visited.front();
		visited.pop();
		visitedSet.insert(i);
	}

	auto exNet = std::make_shared<DNetwork<>>();
	for (auto it = visitedSet.begin(); it != visitedSet.end(); ++it) {
		exNet->addNode(net->getLabel(*it));
	}
	for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it) {
		auto i = net->start(*it);
		auto j = net->end(*it);
		if ((visitedSet.count(i) > 0) && (visitedSet.count(j) > 0)) {
			exNet->addEdge(exNet->getID(net->getLabel(i)),
					exNet->getID(net->getLabel(j)));
		}
	}
	exNet->assemble();
	std::cout << "# Extraced from " << argv[1] << " with srcNode = " << srcNode
			<< " and limit = " << limit << std::endl;
	std::cout << "# Number of nodes: " << exNet->getNbNodes()
			<< ", number of edges: " << exNet->getNbEdges() << std::endl;
	exNet->print();
	return 0;
}

