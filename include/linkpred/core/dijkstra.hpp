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
 * @brief Contains an implementation of Dijkstra's algorithm.
 */

#ifndef DIJKSTRA_HPP_
#define DIJKSTRA_HPP_

#include "LinkPredConfig.hpp"
#include "linkpred/core/unetwork.hpp"
#include "linkpred/core/dnetwork.hpp"
#include <map>
#include <set>
#include <vector>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

namespace LinkPred {

/**
 * @brief An implementation of Dijkstra's algorithm.
 * @details This class currently implements Dijkstra's algorithm using a binary heap,
 * which means that the algorithm runs in O(m + n (log n)^2) instead of O(m + n log n).
 * @tparam NetworkT The network type.
 * @tparam DistType The distance type (can be an integer or floating point type).
 * @tparam NbHopsType The type of the number of hops (must be an integer type).
 */
template<typename NetworkT = UNetwork<>, typename DistType = double,
		typename NbHopsType = std::size_t> class Dijkstra {

public:

	using LengthMapIdType = long int; /**< Length map ID type. */
	using NetworkSP = std::shared_ptr<NetworkT>; /**< Shared pointer to network. */
	using NodeIdType = typename NetworkT::NodeIdType; /**< Nodes ID type. */
	using LabelType = typename NetworkT::LabelType; /**< Nodes label type. */
	using EdgeType = typename NetworkT::EdgeType; /**< Edge type. */
	using NodeDistMap = typename NetworkT::template NodeMap<std::pair<DistType, NbHopsType>>; /**< Distance map. */
	using NodeDistMapSP = typename NetworkT::template NodeMapSP<std::pair<DistType, NbHopsType>>; /**< Shared pointer to a distance map. */
	using NodeSDistMapSP = typename NetworkT::template NodeSMapSP<std::pair<DistType, NbHopsType>>; /**< Shared pointer to a sparse distance map. */
	using EdgeLengthMap = typename NetworkT::template EdgeMap<DistType>; /**< Edge length map. */
	using EdgeLengthMapSP = typename NetworkT::template EdgeMapSP<DistType>; /**< Shared pointer to an edge length map. */
	using PathType = std::vector<NodeIdType>; /**< Path. */
	using PathTypeSP = std::shared_ptr<PathType>; /**< Shared pointer to a path. */

protected:

	std::shared_ptr<NetworkT const> net; /**< The network. */
	std::map<LengthMapIdType, EdgeLengthMapSP> lengthMaps; /**< Length maps of the edges. */
#ifdef WITH_OPENMP
	bool parallel = false; /**< Enable/disable parallelism. */
#endif

	/**
	 * @return a random length map ID.
	 */
	static LengthMapIdType getRndLengthMapId() {
		RandomGen rng;
		return rng.getInt();
	}

public:
	/**
	 * Constructor.
	 */
	Dijkstra(std::shared_ptr<NetworkT const> net) :
			net(net) {
	}

	/**
	 * Default constructor.
	 */
	Dijkstra() = default;

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	Dijkstra(Dijkstra const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	Dijkstra & operator =(Dijkstra const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	Dijkstra(Dijkstra && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	Dijkstra & operator =(Dijkstra && that) = default;

	/**
	 * @return The network.
	 */
	std::shared_ptr<NetworkT const> getNet() const {
		return net;
	}

	/**
	 * Update the network. Alllength maps are invalidated.
	 * @param net
	 */
	void setNet(std::shared_ptr<const NetworkT> net) {
		this->net = net;
		lengthMaps.clear();
	}

#ifdef WITH_OPENMP
	/**
	 * @return Whether parallelism is enabled.
	 */
	bool isParallel() const {
		return parallel;
	}

	/**
	 * Enable/disable parallelism.
	 * @param parallel True to enable parallelism, false to disable it.
	 */
	void setParallel(bool parallel) {
		this->parallel = parallel;
	}
#endif

	/**
	 * Register a length map.
	 * @param length The length map.
	 * @return The ID of the length map.
	 */
	LengthMapIdType registerLengthMap(EdgeLengthMapSP length);

	/**
	 * Unregister a length map.
	 * @param lengthMapId The ID of the length map.
	 */
	void unregisterLengthMap(LengthMapIdType const & lengthMapId);

	/**
	 * Computes the distance between two nodes while ignoring the edge between them.
	 * @param srcId The ID of the source node.
	 * @param dstId The ID of the destination node.
	 * @param lengthMapId The ID of the map length (as returned by UNetwork::registerLengthMap).
	 * @param discDist The value that should be assigned as distance between disconnected nodes.
	 * @param discNbHops The value that should be assigned as number of hops between disconnected nodes.
	 * @return The distance map.
	 */
	std::pair<DistType, NbHopsType> getIndDist(NodeIdType const & srcId,
			NodeIdType const & dstId, LengthMapIdType lengthMapId,
			DistType discDist = std::numeric_limits<DistType>::infinity(),
			NbHopsType discNbHops =
					std::numeric_limits<NbHopsType>::max()) const;

	/**
	 * Computes the distance from a single source node to all other nodes while ignoring the edge between the source node and a specific node.
	 * @param srcId The ID of the source node.
	 * @param dstId The ID of the destination node of the edge to be ignored.
	 * @param lengthMapId The ID of the map length (as returned by UNetwork::registerLengthMap).
	 * @param discDist The value that should be assigned as distance between disconnected nodes.
	 * @param discNbHops The value that should be assigned as number of hops between disconnected nodes.
	 * @return The distance map.
	 */
	NodeDistMapSP getIndDistToAll(NodeIdType const & srcId,
			NodeIdType const & dstId, LengthMapIdType lengthMapId,
			DistType discDist = std::numeric_limits<DistType>::infinity(),
			NbHopsType discNbHops =
					std::numeric_limits<NbHopsType>::max()) const;

	/**
	 * Finds the shortest path from a node to another.
	 * @param srcId The ID of the source node.
	 * @param dstId The ID of the destination node.
	 * @param lengthMapId The ID of the map length (as returned by UNetwork::registerLengthMap).
	 * @param discDist The value that should be assigned as distance between disconnected nodes.
	 * @param discNbHops The value that should be assigned as number of hops between disconnected nodes.
	 * @return A pair containing a pointer to the path and its length.
	 */
	std::pair<PathTypeSP, double> getShortestPath(NodeIdType const & srcId,
			NodeIdType const & dstId, LengthMapIdType lengthMapId,
			DistType discDist = std::numeric_limits<DistType>::infinity(),
			NbHopsType discNbHops =
					std::numeric_limits<NbHopsType>::max()) const;

	/**
	 * Computes the distance from a single source node to all other nodes.
	 * @param srcId The ID of the source node.
	 * @param lengthMapId The ID of the map length (as returned by UNetwork::registerLengthMap).
	 * @param discDist The value that should be assigned as distance between disconnected nodes.
	 * @param discNbHops The value that should be assigned as number of hops between disconnected nodes.
	 * @return The distance map.
	 */
	NodeDistMapSP getDist(NodeIdType const & srcId, LengthMapIdType lengthMapId,
			DistType discDist = std::numeric_limits<DistType>::infinity(),
			NbHopsType discNbHops =
					std::numeric_limits<NbHopsType>::max()) const;

	/**
	 * Computes the distance with a limit on the number of hops  from a single source node to all other nodes.
	 * @param srcId The ID of the source node.
	 * @param lengthMapId The ID of the map length (as returned by UNetwork::registerLengthMap).
	 * @param lim The limit on the number of hops.
	 * @param discDist The value that should be assigned as distance between disconnected nodes.
	 * @param discNbHops The value that should be assigned as number of hops between disconnected nodes.
	 * @return The distance map.
	 */
	NodeSDistMapSP getDistL(NodeIdType const & srcId,
			LengthMapIdType lengthMapId, std::size_t lim, DistType discDist =
					std::numeric_limits<DistType>::infinity(),
			NbHopsType discNbHops =
					std::numeric_limits<NbHopsType>::max()) const;

	/**
	 * Computes the distance from a list of nodes to all other nodes.
	 * @param begin Iterator to the first element in the list of source nodes.
	 * @param end Iterator to one past the last element in the list of source nodes.
	 * @param outit Output iterator.
	 * @param lengthMapId The ID of the map length (as returned by UNetwork::registerLengthMap).
	 * @param discDist The value that should be assigned as distance between disconnected nodes.
	 * @param discNbHops The value that should be assigned as number of hops between disconnected nodes.
	 */
	template<typename InputIterator, typename OuputIterator> void getDist(
			InputIterator begin, InputIterator end, OuputIterator outit,
			LengthMapIdType lengthMapId, DistType discDist =
					std::numeric_limits<DistType>::infinity(),
			NbHopsType discNbHops =
					std::numeric_limits<NbHopsType>::max()) const {

		std::vector<NodeIdType> srcIdsVec;
		for (auto it = begin; it != end; ++it) {
			srcIdsVec.push_back(*it);
		}
		std::vector<NodeDistMapSP> distVec;
		distVec.resize(srcIdsVec.size());
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < srcIdsVec.size(); i++) {
			distVec[i] = getDist(srcIdsVec[i], lengthMapId, discDist,
					discNbHops);
		}

		std::copy(distVec.begin(), distVec.end(), outit);
	}

	/**
	 * Computes the inverse similarity from a list of nodes to all other nodes.
	 * @param begin Iterator to the first element in the list of source nodes.
	 * @param end Iterator to one past the last element in the list of source nodes.
	 * @param outit Output iterator.
	 * @param lengthMapId The ID of the map length (as returned by UNetwork::registerLengthMap).
	 * @param discDist The value that should be assigned as distance between disconnected nodes.
	 * @param discNbHops The value that should be assigned as number of hops between disconnected nodes.
	 */
	template<typename InputIterator, typename OuputIterator> void getSDist(
			InputIterator begin, InputIterator end, OuputIterator outit,
			LengthMapIdType lengthMapId, DistType discDist =
					std::numeric_limits<DistType>::infinity(),
			NbHopsType discNbHops =
					std::numeric_limits<NbHopsType>::max()) const {

		std::vector<NodeIdType> srcIdsVec;
		for (auto it = begin; it != end; ++it) {
			srcIdsVec.push_back(*it);
		}
		std::vector<NodeDistMapSP> distVec;
		distVec.resize(srcIdsVec.size());
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < srcIdsVec.size(); i++) {
			distVec[i] = getSDist(srcIdsVec[i], lengthMapId, discDist,
					discNbHops);
		}

		std::copy(distVec.begin(), distVec.end(), outit);
	}

	/**
	 * Computes the dissimilarity of a single source node to all other nodes.
	 * @param srcId The ID of the source node.
	 * @param lengthMapId The ID of the map length (as returned by UNetwork::registerLengthMap).
	 * @param discDist The value that should be assigned as distance between disconnected nodes.
	 * @param discNbHops The value that should be assigned as number of hops between disconnected nodes.
	 * @return The distance map.
	 */
	NodeDistMapSP getDsim(NodeIdType const & srcId, LengthMapIdType lengthMapId,
			DistType discDist = std::numeric_limits<DistType>::infinity(),
			NbHopsType discNbHops =
					std::numeric_limits<NbHopsType>::max()) const;

	/**
	 * Computes the dissimilarity between two nodes while ignoring the direction edge between them.
	 * @param srcId The ID of the source node.
	 * @param dstId The ID of the destination node.
	 * @param lengthMapId The ID of the map length (as returned by UNetwork::registerLengthMap).
	 * @param discDist The value that should be assigned as distance between disconnected nodes.
	 * @param discNbHops The value that should be assigned as number of hops between disconnected nodes.
	 * @return The distance map.
	 */
	std::pair<DistType, NbHopsType> getIndDsim(NodeIdType const & srcId,
			NodeIdType const & dstId, LengthMapIdType lengthMapId,
			DistType discDist = std::numeric_limits<DistType>::infinity(),
			NbHopsType discNbHops =
					std::numeric_limits<NbHopsType>::max()) const;

	/**
	 * Computes the dissimilarity of a set of nodes to all other nodes.
	 * @param srcIds The IDs of the source node.
	 * @param lengthMapId The ID of the map length (as returned by UNetwork::registerLengthMap).
	 * @param discDist The value that should be assigned as distance between disconnected nodes.
	 * @param discNbHops The value that should be assigned as number of hops between disconnected nodes.
	 * @return The distance map.
	 */
	std::map<NodeIdType, NodeDistMapSP> getDsim(
			std::set<NodeIdType> const & srcIds, LengthMapIdType lengthMapId,
			DistType discDist = std::numeric_limits<DistType>::infinity(),
			NbHopsType discNbHops =
					std::numeric_limits<NbHopsType>::max()) const;

	/**
	 * Destructor.
	 */
	virtual ~Dijkstra() = default;
};

}
/* namespace LinkPred */

#endif /* DIJKSTRA_HPP_ */
