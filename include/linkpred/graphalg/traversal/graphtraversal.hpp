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

/**
 * \file
 * @ingroup GraphAlg
 * @brief Contains the implementation of graph traversal algorithms.
 */

#ifndef GRAPHTRAVERSAL_HPP_
#define GRAPHTRAVERSAL_HPP_

#include "linkpred/core/unetwork/unetwork.hpp"
#include "linkpred/core/dnetwork/dnetwork.hpp"
#include "linkpred/utils/log.hpp"
#include <queue>
#include <stack>
#include <set>

namespace LinkPred {

/**
 * @brief A class that counts nodes during traversal.
 * @tparam Network The network type.
 */
template<typename Network = UNetwork<>> class Counter {
protected:
	std::size_t count = 0; /**< To count the number of nodes visited. */

public:

	/**
	 * Constructor.
	 */
	Counter() = default;

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	Counter(Counter const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	Counter & operator =(Counter const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	Counter(Counter && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	Counter & operator =(Counter && that) = default;

	/**
	 * Node processing.
	 * @param i The node's ID.
	 */
	bool process(typename Network::NodeID const & i) {
		count++;
		return true;
	}

	/**
	 * @return The nodes count.
	 */
	std::size_t getCount() const {
		return count;
	}

	/**
	 * Reset he nodes count to 0.
	 */
	void resetCount() {
		count = 0;
	}

	/**
	 * Destructor.
	 */
	virtual ~Counter() = default;

};

/**
 * @brief A class that collects nodes during traversal.
 * @tparam Network The network type.
 */
template<typename Network = UNetwork<>> class Collector {
protected:
	std::queue<typename Network::NodeID> visited; /**< To store visited nodes. */

public:

	/**
	 * Constructor.
	 */
	Collector() = default;

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	Collector(Collector const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	Collector & operator =(Collector const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	Collector(Collector && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	Collector & operator =(Collector && that) = default;

	/**
	 * Node processing.
	 * @param i The node's ID.
	 */
	bool process(typename Network::NodeID const & i) {
		visited.push(i);
		return true;
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
	virtual ~Collector() = default;

};

/**
 * @brief Graph traversal interface.
 * @tparam Network The network type.
 * @tparam NodeProcessor The node processor type.
 */
template<typename Network = UNetwork<>, typename NodeProcessor = Collector<
		Network>> class GraphTraversal {

protected:
	std::shared_ptr<Network const> net; /**< The network on which traversal is done. */

public:
	/**
	 * Constructor.
	 * @param net The network on which traversal is done.
	 */
	GraphTraversal(std::shared_ptr<Network const> net) :
			net(net) {

	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	GraphTraversal(GraphTraversal const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	GraphTraversal & operator =(GraphTraversal const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	GraphTraversal(GraphTraversal && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	GraphTraversal & operator =(GraphTraversal && that) = default;

	/**
	 * Traverse the graph.
	 * @param srcNode The source node.
	 * @param processor The node processor.
	 */
	virtual void traverse(typename Network::NodeID srcNode,
			NodeProcessor & processor) = 0;

	/**
	 * Destructor.
	 */
	virtual ~GraphTraversal() = default;
};

/**
 * @brief BFS graph traversal.
 * @tparam Network The network type.
 * @tparam NodeProcessor The node processor type.
 */
template<typename Network = UNetwork<>, typename NodeProcessor = Collector<
		Network>> class BFS: public GraphTraversal<Network, NodeProcessor> {

	using GraphTraversal<Network, NodeProcessor>::net; /**< The network on which traversal is done. */

public:
	/**
	 * Constructor.
	 * @param net The network on which traversal is done.
	 */
	BFS(std::shared_ptr<Network const> net) :
			GraphTraversal<Network, NodeProcessor>(net) {

	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	BFS(BFS const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	BFS & operator =(BFS const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	BFS(BFS && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	BFS & operator =(BFS && that) = default;

	/**
	 * Traverse the graph in BFS order.
	 * @param srcNode The source node.
	 * @param processor The node processor.
	 */
	virtual void traverse(typename Network::NodeID srcNode,
			NodeProcessor & processor) {

		std::queue<typename Network::NodeID> q;
		std::set<typename Network::NodeID> visited;
		q.push(srcNode);
		visited.insert(srcNode);

		while (!q.empty()) {
			auto i = q.front();
			q.pop();
			if (!processor.process(i)) {
				break;
			}

			for (auto it = net->neighbBegin(i); it != net->neighbEnd(i); ++it) {
				auto j = Network::end(*it);
				if (visited.count(j) == 0) {
					q.push(j);
					visited.insert(j);
				}
			}
		}
	}

	/**
	 * Destructor.
	 */
	virtual ~BFS() = default;
};

/**
 * @brief DFS graph traversal.
 * @tparam Network The network type.
 * @tparam NodeProcessor The node processor type.
 */
template<typename Network = UNetwork<>, typename NodeProcessor = Collector<
		Network>> class DFS: public GraphTraversal<Network, NodeProcessor> {

	using GraphTraversal<Network, NodeProcessor>::net; /**< The network on which traversal is done. */

public:
	/**
	 * Constructor.
	 * @param net The network on which traversal is done.
	 */
	DFS(std::shared_ptr<Network const> net) :
			GraphTraversal<Network, NodeProcessor>(net) {

	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	DFS(DFS const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	DFS & operator =(DFS const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	DFS(DFS && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	DFS & operator =(DFS && that) = default;

	/**
	 * Traverse the graph in DFS order.
	 * @param srcNode The source node.
	 * @param processor The node processor.
	 */
	virtual void traverse(typename Network::NodeID srcNode,
			NodeProcessor & processor) {

		std::stack<typename Network::NodeID> st;
		std::set<typename Network::NodeID> visited;
		st.push(srcNode);
		visited.insert(srcNode);

		while (!st.empty()) {
			auto i = st.top();
			st.pop();
			if (!processor.process(i)) {
				break;
			}

			for (auto it = net->neighbBegin(i); it != net->neighbEnd(i); ++it) {
				auto j = Network::end(*it);
				if (visited.count(j) == 0) {
					st.push(j);
					visited.insert(j);
				}
			}
		}
	}

	/**
	 * Destructor.
	 */
	virtual ~DFS() = default;
};

} /* namespace LinkPred */

#endif /* INCLUDE_GRAPHTRAVERSAL_HPP_ */
