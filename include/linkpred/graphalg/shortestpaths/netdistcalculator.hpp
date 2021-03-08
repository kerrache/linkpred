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
 * @brief Contains the implementation of classes for computing distances in network.
 */

#ifndef DISTCALCULATOR_HPP_
#define DISTCALCULATOR_HPP_

#include "linkpred/core/unetwork/unetwork.hpp"
#include "linkpred/core/dnetwork/dnetwork.hpp"
#include "linkpred/graphalg/shortestpaths/dijkstra.hpp"
#include <tuple>
#include <queue>

namespace LinkPred {

/**
 * @brief Cache levels.
 */
enum CacheLevel {
	NoCache, /**< No cache. */
	NodeCache, /**< Cache at node level. */
	NetworkCache /**< Cache at network level. */
};

/**
 * @brief Interface for calculating the distance between nodes in a network.
 * @tparam Network The network type.
 * @tparam DistType The distance type (can be an integer or floating point type).
 * @tparam NbHopsType The type of the number of hops (must be an integer type).
 */
template<typename Network, typename DistType, typename NbHopsType> class NetDistCalculator {

public:

	using LengthMapIdType = long int; /**< Length map ID type. */
	using NetworkSP = std::shared_ptr<Network>; /**< Shared pointer to network. */
	using NodeID = typename Network::NodeID; /**< Nodes ID type. */
	using Label = typename Network::Label; /**< Nodes label type. */
	using Edge = typename Network::Edge; /**< Edge type. */
	using NodeDistMap = typename Network::template NodeMap<std::pair<DistType, NbHopsType>>; /**< Distance map. */
	using NodeDistMapSP = typename Network::template NodeMapSP<std::pair<DistType, NbHopsType>>; /**< Shared pointer to a distance map. */
	using EdgeLengthMap = typename Network::template EdgeMap<DistType>; /**< Edge length map. */
	using EdgeLengthMapSP = typename Network::template EdgeMapSP<DistType>; /**< Shared pointer to an edge length map. */

protected:
	DistType discDist = std::numeric_limits<DistType>::infinity(); /**< The value that should be assigned as distance between disconnected nodes. */
	NbHopsType discNbHops = std::numeric_limits<NbHopsType>::max(); /**< The value that should be assigned as number of hops between disconnected nodes. */

public:

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The distance between i and j.
	 */
	virtual std::pair<DistType, NbHopsType> getDist(NodeID const & i,
			NodeID const & j) = 0;

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The distance between i and j ignoring the edge between i and j.
	 */
	virtual std::pair<DistType, NbHopsType> getIndDist(NodeID const & i,
			NodeID const & j) = 0;

	/**
	 * @param i Source node.
	 * @return The distance from i to all other nodes.
	 */
	virtual NodeDistMapSP getDist(NodeID const & i) = 0;

	/**
	 * @param i Index of the source node.
	 * @return The maximum distance from i to any other node in the same connected component.
	 */
	virtual std::pair<DistType, NbHopsType> getMaxDist(NodeID const & i) {
//		std::cerr << "MaxDist..." << std::endl;
		auto maxDist = std::make_pair(discDist, discNbHops);
		auto dm = this->getDist(i);
		for (auto it = dm->begin(); it != dm->end(); ++it) {
			auto dist = *it;
			if (dist.first < maxDist.first) {
				maxDist = dist;
			}
		}
		if (maxDist.second == discNbHops) {
			maxDist = std::make_pair(0, 0);
		}
//		std::cerr << "Done..." << std::endl;
		return maxDist;
	}

	/**
	 * @return The distance assigned to disconnected nodes.
	 */
	double getDiscDist() const {
		return discDist;
	}

	/**
	 * Set the distance assigned to disconnected nodes.
	 * @param discDist The distance assigned to disconnected nodes.
	 */
	void setDiscDist(DistType discDist) {
		this->discDist = discDist;
	}

	/**
	 * @return The value that should be assigned as number of hops between disconnected nodes.
	 */
	std::size_t getDiscNbHops() const {
		return discNbHops;
	}

	/**
	 * Set the value that should be assigned as number of hops between disconnected nodes.
	 * @param discNbHops The value that should be assigned as number of hops between disconnected nodes.
	 */
	void setDiscNbHops(NbHopsType discNbHops) {
		this->discNbHops = discNbHops;
	}

	/**
	 * Destructor.
	 */
	virtual ~NetDistCalculator() = default;
};

/**
 * @brief Interface for calculating the similarity between nodes in a network.
 * @tparam Network The network type.
 * @tparam DistType The distance and similarity type (can be an integer or floating point type).
 * @tparam NbHopsType The type of the number of hops (must be an integer type).
 */
template<typename Network, typename DistType, typename NbHopsType> class NetSimlCalculator {

public:

	using LengthMapIdType = long int; /**< Length map ID type. */
	using NetworkSP = std::shared_ptr<Network>; /**< Shared pointer to network. */
	using NodeID = typename Network::NodeID; /**< Nodes ID type. */
	using Label = typename Network::Label; /**< Nodes label type. */
	using Edge = typename Network::Edge; /**< Edge type. */
	using EdgeLengthMap = typename Network::template EdgeMap<DistType>; /**< Edge length map. */
	using EdgeLengthMapSP = typename Network::template EdgeMapSP<DistType>; /**< Shared pointer to an edge length map. */
	using NodeSDistMap = typename Network::template NodeSMap<std::pair<DistType, NbHopsType>>; /**< Distance map. */
	using NodeSDistMapSP = typename Network::template NodeSMapSP<std::pair<DistType, NbHopsType>>; /**< Shared pointer to a distance map. */

protected:
	DistType selfSiml = std::numeric_limits<DistType>::infinity(); /**< The value that should be assigned as similarity between a node and itself. */
	NbHopsType discNbHops = std::numeric_limits<NbHopsType>::max(); /**< The value that should be assigned as number of hops between disconnected nodes. */
	DistType discDist = std::numeric_limits<DistType>::infinity(); /**< The value that should be assigned as distance between disconnected nodes. */
	double lambda = 0.5; /**< Conductance. */

public:

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The similarity between i and j.
	 */
	virtual std::tuple<DistType, DistType, NbHopsType> getSiml(
			NodeID const & i, NodeID const & j) = 0;

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The directed similarity between i and j.
	 */
	virtual std::tuple<DistType, DistType, NbHopsType> getDirSiml(
			NodeID const & i, NodeID const & j) = 0;

	/**
	 * @return A sparse similarity map of nodes having non-zero similarity to a given source node. Only nodes not connected to srcNode ar considered.
	 * @param srcNode The source node.
	 */
	virtual NodeSDistMapSP getNnzSimlMap(NodeID const & srcNode) {
		throw std::runtime_error("Not implemented");
	}

	/**
	 * @return A sparse similarity map of nodes having non-zero similarity to a given source node.
	 * Only nodes not connected to srcNode are considered. The node srcNode itself is also excluded.
	 * @param srcNode The source node.
	 */
	virtual NodeSDistMapSP getNnzSimlMapNoNeighb(NodeID const & srcNode) {
		throw std::runtime_error("Not implemented");
	}

	/**
	 * @return The similarity between a node and itself.
	 */
	double getSelfSiml() const {
		return selfSiml;
	}

	/**
	 * Set the similarity between a node and itself.
	 * @param selfSiml The similarity similarity between a node and itself.
	 */
	void setSelfSiml(DistType selfSiml) {
		this->selfSiml = selfSiml;
	}

	/**
	 * @return The distance assigned to disconnected nodes.
	 */
	double getDiscDist() const {
		return discDist;
	}

	/**
	 * Set the distance assigned to disconnected nodes.
	 * @param discDist The distance assigned to disconnected nodes.
	 */
	void setDiscDist(DistType discDist) {
		this->discDist = discDist;
	}

	/**
	 * @return The value that should be assigned as number of hops between disconnected nodes.
	 */
	std::size_t getDiscNbHops() const {
		return discNbHops;
	}

	/**
	 * Set the value that should be assigned as number of hops between disconnected nodes.
	 * @param discNbHops The value that should be assigned as number of hops between disconnected nodes.
	 */
	void setDiscNbHops(NbHopsType discNbHops) {
		this->discNbHops = discNbHops;
	}

	/**
	 * @return The conductance.
	 */
	double getLambda() const {
		return lambda;
	}

	/**
	 * Set the conductance.
	 * @param lambda The new value of the conductance.
	 */
	void setLambda(double lambda) {
		this->lambda = lambda;
	}

	/**
	 * Destructor.
	 */
	virtual ~NetSimlCalculator() = default;
};

/**
 * @brief Interface for calculating the indirect similarity between nodes in a network.
 * @tparam Network The network type.
 * @tparam DistType The distance and similarity type (can be an integer or floating point type).
 * @tparam NbHopsType The type of the number of hops (must be an integer type).
 */
template<typename Network, typename DistType, typename NbHopsType> class NetIndSimlCalculator {

public:

	using LengthMapIdType = long int; /**< Length map ID type. */
	using NetworkSP = std::shared_ptr<Network>; /**< Shared pointer to network. */
	using NodeID = typename Network::NodeID; /**< Nodes ID type. */
	using Label = typename Network::Label; /**< Nodes label type. */
	using Edge = typename Network::Edge; /**< Edge type. */
	using EdgeLengthMap = typename Network::template EdgeMap<DistType>; /**< Edge length map. */
	using EdgeLengthMapSP = typename Network::template EdgeMapSP<DistType>; /**< Shared pointer to an edge length map. */

protected:
	DistType selfSiml = std::numeric_limits<DistType>::infinity(); /**< The value that should be assigned as similarity between a node and itself. */
	NbHopsType discNbHops = std::numeric_limits<NbHopsType>::max(); /**< The value that should be assigned as number of hops between disconnected nodes. */
	DistType discDist = std::numeric_limits<DistType>::infinity(); /**< The value that should be assigned as distance between disconnected nodes. */
	double lambda = 0.5; /**< Balance between distance and number of hops. */

public:

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The similarity between i and j ignoring the edge between i and j.
	 */
	virtual std::tuple<DistType, DistType, NbHopsType> getIndSiml(
			NodeID const & i, NodeID const & j) = 0;

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The directed similarity between i and j ignoring the edge between i and j.
	 */
	virtual std::tuple<DistType, DistType, NbHopsType> getDirIndSiml(
			NodeID const & i, NodeID const & j) = 0;

	/**
	 * @return The similarity between a node and itself.
	 */
	double getSelfSiml() const {
		return selfSiml;
	}

	/**
	 * Set the similarity between a node and itself.
	 * @param selfSiml The similarity similarity between a node and itself.
	 */
	void setSelfSiml(DistType selfSiml) {
		this->selfSiml = selfSiml;
	}

	/**
	 * @return The distance assigned to disconnected nodes.
	 */
	double getDiscDist() const {
		return discDist;
	}

	/**
	 * Set the distance assigned to disconnected nodes.
	 * @param discDist The distance assigned to disconnected nodes.
	 */
	void setDiscDist(DistType discDist) {
		this->discDist = discDist;
	}

	/**
	 * @return The value that should be assigned as number of hops between disconnected nodes.
	 */
	std::size_t getDiscNbHops() const {
		return discNbHops;
	}

	/**
	 * Set the value that should be assigned as number of hops between disconnected nodes.
	 * @param discNbHops The value that should be assigned as number of hops between disconnected nodes.
	 */
	void setDiscNbHops(NbHopsType discNbHops) {
		this->discNbHops = discNbHops;
	}

	/**
	 * @return The conductance.
	 */
	double getLambda() const {
		return lambda;
	}

	/**
	 * Set the conductance.
	 * @param lambda The new value of the weight.
	 */
	void setLambda(double lambda) {
		this->lambda = lambda;
	}

	/**
	 * Destructor.
	 */
	virtual ~NetIndSimlCalculator() = default;
};

/**
 * @brief Exact shortest path distance calculator.
 * @details This class offers an additional layer over Dijkstra which provides
 * memory management functionalities by caching computed distances.
 * @tparam Network The network type.
 * @tparam DistType The distance type (can be an integer or floating point type).
 * @tparam NbHopsType The type of the number of hops (must be an integer type).
 */
template<typename Network = UNetwork<>, typename DistType = double,
		typename NbHopsType = std::size_t> class ESPDistCalculator: public NetDistCalculator<
		Network, DistType, NbHopsType> {

public:
	using LengthMapIdType= typename NetDistCalculator<Network, DistType, NbHopsType>::LengthMapIdType; /**< Length map ID type. */
	using NetworkSP = typename NetDistCalculator<Network, DistType, NbHopsType>::NetworkSP; /**< Shared pointer to network. */
	using NodeID = typename NetDistCalculator<Network, DistType, NbHopsType>::NodeID; /**< Nodes ID type. */
	using Label = typename NetDistCalculator<Network, DistType, NbHopsType>::Label; /**< Nodes label type. */
	using Edge = typename NetDistCalculator<Network, DistType, NbHopsType>::Edge; /**< Edge type. */
	using NodeDistMap = typename NetDistCalculator<Network, DistType, NbHopsType>::NodeDistMap; /**< Distance map. */
	using NodeDistMapSP = typename NetDistCalculator<Network, DistType, NbHopsType>::NodeDistMapSP; /**< Shared pointer to a distance map. */
	using EdgeLengthMap = typename NetDistCalculator<Network, DistType, NbHopsType>::EdgeLengthMap; /**< Edge length map. */
	using EdgeLengthMapSP = typename NetDistCalculator<Network, DistType, NbHopsType>::EdgeLengthMapSP; /**< Shared pointer to an edge length map. */

protected:
	Dijkstra<Network, DistType, NbHopsType> & dijkstra; /**< Dijkstra's algorithm. */
	EdgeLengthMapSP length; /**< The length map. */
	LengthMapIdType lengthMapId; /**< The ID of the length map. */
	CacheLevel cacheLevel; /**< The level of distance cache. */
	NodeDistMapSP nodeDistCache; /**< Cache of distances from a node to all other nodes.*/
	NodeID cachedNodeId; /**< Previous node. */
	std::map<NodeID, NodeDistMapSP> netDistCache; /**< Cache of distances from a node to all other nodes.*/
	std::size_t cacheHits = 0; /**< Number of cache hits. */
	std::size_t cacheMiss = 0; /**< Number of cahce misses. */

public:
	/**
	 * Constructor.
	 * @param dijkstra A Dijkstra algorithm object.
	 * @param length The length map.
	 * @param cacheLevel The cache level.
	 */
	ESPDistCalculator(Dijkstra<Network, DistType, NbHopsType> & dijkstra,
			EdgeLengthMapSP length, CacheLevel cacheLevel =
					CacheLevel::NetworkCache) :
			dijkstra(dijkstra), length(length), cacheLevel(cacheLevel) {

		lengthMapId = dijkstra.registerLengthMap(length);
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	ESPDistCalculator(ESPDistCalculator const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	ESPDistCalculator & operator =(ESPDistCalculator const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	ESPDistCalculator(ESPDistCalculator && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	ESPDistCalculator & operator =(ESPDistCalculator && that) = default;

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The distance between i and j.
	 */
	virtual std::pair<DistType, NbHopsType> getDist(NodeID const & i,
			NodeID const & j);

	/**
	 * @param i Source node.
	 * @return The distance from i to all other nodes.
	 */
	virtual NodeDistMapSP getDist(NodeID const & i);

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The distance between i and j ignoring the edge between i and j.
	 */
	virtual std::pair<DistType, NbHopsType> getIndDist(NodeID const & i,
			NodeID const & j);

	/**
	 * @return The number of cache hits.
	 */
	std::size_t getCacheHits() const {
		return cacheHits;
	}

	/**
	 * @return The number of cache misses.
	 */
	std::size_t getCacheMiss() const {
		return cacheMiss;
	}

	/**
	 * Destructor.
	 */
	virtual ~ESPDistCalculator() {
		dijkstra.unregisterLengthMap(lengthMapId);
	}
};

/**
 * @brief Exact shortest path distance calculator with limits on the number of hops.
 * @details This class offers an additional layer over Dijkstra which provides
 * memory management functionalities by caching computed distances.
 * @tparam Network The network type.
 * @tparam DistType The distance type (can be an integer or floating point type).
 * @tparam NbHopsType The type of the number of hops (must be an integer type).
 */
template<typename Network = UNetwork<>, typename DistType = double,
		typename NbHopsType = std::size_t> class ESPLDistCalculator: public NetDistCalculator<
		Network, DistType, NbHopsType> {

public:
	using LengthMapIdType= typename NetDistCalculator<Network, DistType, NbHopsType>::LengthMapIdType; /**< Length map ID type. */
	using NetworkSP = typename NetDistCalculator<Network, DistType, NbHopsType>::NetworkSP; /**< Shared pointer to network. */
	using NodeID = typename NetDistCalculator<Network, DistType, NbHopsType>::NodeID; /**< Nodes ID type. */
	using Label = typename NetDistCalculator<Network, DistType, NbHopsType>::Label; /**< Nodes label type. */
	using Edge = typename NetDistCalculator<Network, DistType, NbHopsType>::Edge; /**< Edge type. */
	using NodeDistMap = typename NetDistCalculator<Network, DistType, NbHopsType>::NodeDistMap; /**< Distance map. */
	using NodeDistMapSP = typename NetDistCalculator<Network, DistType, NbHopsType>::NodeDistMapSP; /**< Shared pointer to a distance map. */
	using NodeSDistMap = typename Network::template NodeSMap<std::pair<DistType, NbHopsType>>; /**< Distance map. */
	using NodeSDistMapSP = typename Network::template NodeSMapSP<std::pair<DistType, NbHopsType>>; /**< Shared pointer to a distance map. */
	using EdgeLengthMap = typename NetDistCalculator<Network, DistType, NbHopsType>::EdgeLengthMap; /**< Edge length map. */
	using EdgeLengthMapSP = typename NetDistCalculator<Network, DistType, NbHopsType>::EdgeLengthMapSP; /**< Shared pointer to an edge length map. */

protected:
	using NetDistCalculator<Network, DistType, NbHopsType>::discNbHops; /**< The value that should be assigned as number of hops between disconnected nodes. */
	using NetDistCalculator<Network, DistType, NbHopsType>::discDist; /**< The value that should be assigned as distance between disconnected nodes. */
	Dijkstra<Network, DistType, NbHopsType> & dijkstra; /**< Dijkstra's algorithm. */
	EdgeLengthMapSP length; /**< The length map. */
	std::size_t lim; /**< Limit of number of hops. */
	std::shared_ptr<Network const> net; /**< The network. */
	LengthMapIdType lengthMapId; /**< The ID of the length map. */
	CacheLevel cacheLevel; /**< The level of distance cache. */
	NodeSDistMapSP nodeDistCache; /**< Cache of distances from a node to all other nodes.*/
	std::queue<NodeID> cachedNodeIds; /**< Nodes in cache. */
	std::size_t maxNbNodesInCache; /**< Maximum number of nodes kept in cache (by default 10% of the nodes). */
	std::map<NodeID, NodeSDistMapSP> netDistCache; /**< Cache of distance of a node to all other nodes.*/

public:
	/**
	 * Constructor.
	 * @param dijkstra A Dijkstra algorithm object.
	 * @param length The length map.
	 * @param lim Horizon limit.
	 * @param cacheLevel The cache level.
	 */
	ESPLDistCalculator(Dijkstra<Network, DistType, NbHopsType> & dijkstra,
			EdgeLengthMapSP length, std::size_t lim, CacheLevel cacheLevel =
					CacheLevel::NetworkCache) :
			dijkstra(dijkstra), length(length), lim(lim), cacheLevel(cacheLevel) {
		net = dijkstra.getNet();
		lengthMapId = dijkstra.registerLengthMap(length);
		maxNbNodesInCache = net->getNbNodes() / 10;
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	ESPLDistCalculator(ESPLDistCalculator const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	ESPLDistCalculator & operator =(ESPLDistCalculator const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	ESPLDistCalculator(ESPLDistCalculator && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	ESPLDistCalculator & operator =(ESPLDistCalculator && that) = default;

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The distance between i and j.
	 */
	virtual std::pair<DistType, NbHopsType> getDist(NodeID const & i,
			NodeID const & j);

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The distance between i and j ignoring the edge between i and j.
	 */
	virtual std::pair<DistType, NbHopsType> getIndDist(NodeID const & i,
			NodeID const & j) {
		throw std::runtime_error("Not implemented");
	}

	/**
	 * @return The distance from srcNode.
	 * @param i The source node.
	 */
	virtual NodeDistMapSP getDist(NodeID const & i);

	/**
	 * @return The sparse distance from srcNode.
	 * @param i The source node.
	 */
	virtual NodeSDistMapSP getSDistMap(NodeID const & i);

	/**
	 * @return A sparse distance map of nodes having finite distance to a given source node.
	 * Only nodes not connected to srcNode are considered. The node srcNode itself is also excluded.
	 * @param srcNode The source node.
	 */
	virtual NodeSDistMapSP getFinDistMapNoNeighb(NodeID const & srcNode);

	/**
	 * @param e An input edge.
	 * @return The distance between start and end of e.
	 */
	virtual std::pair<DistType, NbHopsType> getDist(Edge const & e) {
		return getDist(net->start(e), net->end(e));
	}

	/**
	 * @return The maximum number of nodes allowed in cache.
	 */
	std::size_t getMaxNbNodesInCache() const {
		return maxNbNodesInCache;
	}

	/**
	 * Set the maximum number of nodes in cache (for each node a map of distance to all other nodes is kept in memory).
	 * @param maxNbNodesInCache New value for the maximum number of nodes in cache.
	 */
	void setMaxNbNodesInCache(std::size_t maxNbNodesInCache) {
		this->maxNbNodesInCache = maxNbNodesInCache;
	}

	/**
	 * @return The limit on the number of hops.
	 */
	std::size_t getLim() const {
		return lim;
	}

	/**
	 * Destructor.
	 */
	virtual ~ ESPLDistCalculator() {
		dijkstra.unregisterLengthMap(lengthMapId);
	}

};

/**
 * @brief Approximate shortest path distance calculator.
 * @details This class implements a fast method to approximate shortest path distances
 * in large networks. First, a set of landmarks is provided, then the distance between
 * any two nodes is computed as the minimum over the sums of distances between the two
 * nodes and any landmark.
 * @tparam Network The network type.
 * @tparam DistType The distance type (can be an integer or floating point type).
 * @tparam NbHopsType The type of the number of hops (must be an integer type).
 */
template<typename Network = UNetwork<>, typename DistType = double,
		typename NbHopsType = std::size_t> class ASPDistCalculator: public NetDistCalculator<
		Network, DistType, NbHopsType> {

public:
	using LengthMapIdType= typename NetDistCalculator<Network, DistType, NbHopsType>::LengthMapIdType; /**< Length map ID type. */
	using NetworkSP = typename NetDistCalculator<Network, DistType, NbHopsType>::NetworkSP; /**< Shared pointer to network. */
	using NodeID = typename NetDistCalculator<Network, DistType, NbHopsType>::NodeID; /**< Nodes ID type. */
	using Label = typename NetDistCalculator<Network, DistType, NbHopsType>::Label; /**< Nodes label type. */
	using Edge = typename NetDistCalculator<Network, DistType, NbHopsType>::Edge; /**< Edge type. */
	using NodeDistMap = typename NetDistCalculator<Network, DistType, NbHopsType>::NodeDistMap; /**< Distance map. */
	using NodeDistMapSP = typename NetDistCalculator<Network, DistType, NbHopsType>::NodeDistMapSP; /**< Shared pointer to a distance map. */
	using EdgeLengthMap = typename NetDistCalculator<Network, DistType, NbHopsType>::EdgeLengthMap; /**< Edge length map. */
	using EdgeLengthMapSP = typename NetDistCalculator<Network, DistType, NbHopsType>::EdgeLengthMapSP; /**< Shared pointer to an edge length map. */

protected:
	Dijkstra<Network, DistType, NbHopsType> & dijkstra; /**< Dijkstra's algorithm. */
	EdgeLengthMapSP length; /**< The length map. */
	LengthMapIdType lengthMapId; /**< The ID of the length map. */
	std::map<NodeID, NodeDistMapSP> distToLandmarks; /**< Distances to landmarks. */

public:
	/**
	 * Constructor.
	 * @param dijkstra A Dijkstra algorithm object.
	 * @param length The length map.
	 */
	ASPDistCalculator(Dijkstra<Network, DistType, NbHopsType> & dijkstra,
			EdgeLengthMapSP length) :
			dijkstra(dijkstra), length(length) {

		lengthMapId = dijkstra.registerLengthMap(length);
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	ASPDistCalculator(ASPDistCalculator const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	ASPDistCalculator & operator =(ASPDistCalculator const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	ASPDistCalculator(ASPDistCalculator && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	ASPDistCalculator & operator =(ASPDistCalculator && that) = default;

	/**
	 * Set the landmarks.
	 * @param landmarksBegin Iterator to the first landmarks.
	 * @param landmarksEnd Iterator to one past the last landmarks.
	 */
	template<typename InputIterator> void setLandmarks(
			InputIterator landmarksBegin, InputIterator landmarksEnd) {

		distToLandmarks.clear();
		std::vector<NodeDistMapSP> dists;
		dijkstra.getDist(landmarksBegin, landmarksEnd,
				std::back_inserter(dists), lengthMapId,
				NetDistCalculator<Network, DistType, NbHopsType>::discDist,
				NetDistCalculator<Network, DistType, NbHopsType>::discNbHops);
		auto dit = dists.begin();
		for (auto it = landmarksBegin; it != landmarksEnd; ++it, ++dit) {
			distToLandmarks[*it] = *dit;
		}
	}

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The distance between i and j.
	 */
	virtual std::pair<DistType, NbHopsType> getDist(NodeID const & i,
			NodeID const & j);

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The distance between i and j ignoring the edge between i and j.
	 */
	virtual std::pair<DistType, NbHopsType> getIndDist(NodeID const & i,
			NodeID const & j) {
		throw std::runtime_error("Not implemented");
	}

	/**
	 * @param i Source node.
	 * @return The distance from i to all other nodes.
	 */
	virtual NodeDistMapSP getDist(NodeID const & i);

	/**
	 * Destructor.
	 */
	virtual ~ASPDistCalculator() {
		dijkstra.unregisterLengthMap(lengthMapId);
	}
};

/**
 * @brief Exact shortest path dissimilarity calculator.
 * @tparam Network The network type.
 * @tparam DistType The distance type (can be an integer or floating point type).
 * @tparam NbHopsType The type of the number of hops (must be an integer type).
 */
template<typename Network = UNetwork<>, typename DsimType = double,
		typename NbHopsType = std::size_t> class ESPDsimCalculator: public NetDistCalculator<
		Network, DsimType, NbHopsType> {

public:
	using LengthMapIdType= typename NetDistCalculator<Network, DsimType, NbHopsType>::LengthMapIdType; /**< Length map ID type. */
	using NetworkSP = typename NetDistCalculator<Network, DsimType, NbHopsType>::NetworkSP; /**< Shared pointer to network. */
	using NodeID = typename NetDistCalculator<Network, DsimType, NbHopsType>::NodeID; /**< Nodes ID type. */
	using Label = typename NetDistCalculator<Network, DsimType, NbHopsType>::Label; /**< Nodes label type. */
	using Edge = typename NetDistCalculator<Network, DsimType, NbHopsType>::Edge; /**< Edge type. */
	using NodeDsimMap = typename NetDistCalculator<Network, DsimType, NbHopsType>::NodeDistMap; /**< Dissimilarity map. */
	using NodeDsimMapSP = typename NetDistCalculator<Network, DsimType, NbHopsType>::NodeDistMapSP; /**< Shared pointer to a Dissimilarity map. */
	using EdgeLengthMap = typename NetDistCalculator<Network, DsimType, NbHopsType>::EdgeLengthMap; /**< Edge length map. */
	using EdgeLengthMapSP = typename NetDistCalculator<Network, DsimType, NbHopsType>::EdgeLengthMapSP; /**< Shared pointer to an edge length map. */

protected:
	Dijkstra<Network, DsimType, NbHopsType> & dijkstra; /**< Dijkstra's algorithm. */
	EdgeLengthMapSP length; /**< The length map. */
	LengthMapIdType lengthMapId; /**< The ID of the length map. */
	CacheLevel cacheLevel; /**< The level of distance cache. */
	NodeDsimMapSP nodeDsimCache; /**< Cache of distances from a node to all other nodes.*/
	NodeID cachedNodeId; /**< Previous node. */
	std::map<NodeID, NodeDsimMapSP> netDsimCache; /**< Cache of distances from a node to all other nodes.*/

public:
	/**
	 * Constructor.
	 * @param dijkstra A Dijkstra algorithm object.
	 * @param length The length map.
	 * @param cacheLevel The cache level.
	 */
	ESPDsimCalculator(Dijkstra<Network, DsimType, NbHopsType> & dijkstra,
			EdgeLengthMapSP length, CacheLevel cacheLevel =
					CacheLevel::NodeCache) :
			dijkstra(dijkstra), length(length), cacheLevel(cacheLevel) {

		lengthMapId = dijkstra.registerLengthMap(length);
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	ESPDsimCalculator(ESPDsimCalculator const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	ESPDsimCalculator & operator =(ESPDsimCalculator const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	ESPDsimCalculator(ESPDsimCalculator && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	ESPDsimCalculator & operator =(ESPDsimCalculator && that) = default;

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The distance between i and j.
	 */
	virtual std::pair<DsimType, NbHopsType> getDist(NodeID const & i,
			NodeID const & j);

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The distance between i and j ignoring the edge between i and j.
	 */
	virtual std::pair<DsimType, NbHopsType> getIndDist(NodeID const & i,
			NodeID const & j);

	/**
	 * @param i Source node.
	 * @return The distance from i to all other nodes.
	 */
	virtual NodeDsimMapSP getDist(NodeID const & i);

	/**
	 * Destructor.
	 */
	virtual ~ESPDsimCalculator() {
		dijkstra.unregisterLengthMap(lengthMapId);
	}
};

/**
 * @brief Approximate shortest path dissimilarity calculator.
 * @tparam Network The network type.
 * @tparam DistType The distance type (can be an integer or floating point type).
 * @tparam NbHopsType The type of the number of hops (must be an integer type).
 */
template<typename Network = UNetwork<>, typename DsimType = double,
		typename NbHopsType = std::size_t> class ASPDsimCalculator: public NetDistCalculator<
		Network, DsimType, NbHopsType> {

public:
	using LengthMapIdType= typename NetDistCalculator<Network, DsimType, NbHopsType>::LengthMapIdType; /**< Length map ID type. */
	using NetworkSP = typename NetDistCalculator<Network, DsimType, NbHopsType>::NetworkSP; /**< Shared pointer to network. */
	using NodeID = typename NetDistCalculator<Network, DsimType, NbHopsType>::NodeID; /**< Nodes ID type. */
	using Label = typename NetDistCalculator<Network, DsimType, NbHopsType>::Label; /**< Nodes label type. */
	using Edge = typename NetDistCalculator<Network, DsimType, NbHopsType>::Edge; /**< Edge type. */
	using NodeDsimMap = typename NetDistCalculator<Network, DsimType, NbHopsType>::NodeDistMap; /**< Dissimilarity map. */
	using NodeDsimMapSP = typename NetDistCalculator<Network, DsimType, NbHopsType>::NodeDistMapSP; /**< Shared pointer to a Dissimilarity map. */
	using EdgeLengthMap = typename NetDistCalculator<Network, DsimType, NbHopsType>::EdgeLengthMap; /**< Edge length map. */
	using EdgeLengthMapSP = typename NetDistCalculator<Network, DsimType, NbHopsType>::EdgeLengthMapSP; /**< Shared pointer to an edge length map. */

protected:
	Dijkstra<Network, DsimType, NbHopsType> & dijkstra; /**< Dijkstra's algorithm. */
	EdgeLengthMapSP length; /**< The length map. */
	LengthMapIdType lengthMapId; /**< The ID of the length map. */
	std::map<NodeID, NodeDsimMapSP> dsimToLandmarks; /**< distances to landmarks. */

public:
	/**
	 * Constructor.
	 * @param dijkstra A Dijkstra algorithm object.
	 * @param length The length map.
	 */
	ASPDsimCalculator(Dijkstra<Network, DsimType, NbHopsType> & dijkstra,
			EdgeLengthMapSP length) :
			dijkstra(dijkstra), length(length) {

		lengthMapId = dijkstra.registerLengthMap(length);
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	ASPDsimCalculator(ASPDsimCalculator const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	ASPDsimCalculator & operator =(ASPDsimCalculator const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	ASPDsimCalculator(ASPDsimCalculator && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	ASPDsimCalculator & operator =(ASPDsimCalculator && that) = default;

	/**
	 * Set the landmarks.
	 * @param landmarksBegin Iterator to the first landmarks.
	 * @param landmarksEnd Iterator to one past the last landmarks.
	 */
	template<typename InputIterator> void setLandmarks(
			InputIterator landmarksBegin, InputIterator landmarksEnd) {
		dsimToLandmarks.clear();
		std::vector<NodeDsimMapSP> dsims;
		dijkstra.getDist(landmarksBegin, landmarksEnd,
				std::back_inserter(dsims), lengthMapId,
				NetDistCalculator<Network, DsimType, NbHopsType>::discDist,
				NetDistCalculator<Network, DsimType, NbHopsType>::discNbHops);
		auto dit = dsims.begin();
		for (auto it = landmarksBegin; it != landmarksEnd; ++it, ++dit) {
			dsimToLandmarks[*it] = *dit;
		}
	}

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The distance between i and j.
	 */
	virtual std::pair<DsimType, NbHopsType> getDist(NodeID const & i,
			NodeID const & j);

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The distance between i and j ignoring the edge between i and j.
	 */
	virtual std::pair<DsimType, NbHopsType> getIndDist(NodeID const & i,
			NodeID const & j) {
		throw std::runtime_error("Not implemented");
	}

	/**
	 * @param i Source node.
	 * @return The distance from i to all other nodes.
	 */
	virtual NodeDsimMapSP getDist(NodeID const & i);

	/**
	 * Destructor.
	 */
	virtual ~ASPDsimCalculator() {
		dijkstra.unregisterLengthMap(lengthMapId);
	}
};

/*******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 *******************************************************************************
 */
/**
 * @brief Exact shortest path distance calculator.
 * @details This class offers an additional layer over Dijkstra which provides
 * memory management functionalities by caching computed distances.
 * @tparam Network The network type.
 * @tparam DistType The distance type (can be an integer or floating point type).
 * @tparam NbHopsType The type of the number of hops (must be an integer type).
 */
template<typename Network = UNetwork<>, typename DistType = double,
		typename NbHopsType = std::size_t> class ESPSimlCalculator: public NetSimlCalculator<
		Network, DistType, NbHopsType> {

public:
	using LengthMapIdType= typename NetSimlCalculator<Network, DistType, NbHopsType>::LengthMapIdType; /**< Length map ID type. */
	using NetworkSP = typename NetSimlCalculator<Network, DistType, NbHopsType>::NetworkSP; /**< Shared pointer to network. */
	using NodeID = typename NetSimlCalculator<Network, DistType, NbHopsType>::NodeID; /**< Nodes ID type. */
	using Label = typename NetSimlCalculator<Network, DistType, NbHopsType>::Label; /**< Nodes label type. */
	using Edge = typename NetSimlCalculator<Network, DistType, NbHopsType>::Edge; /**< Edge type. */
	using NodeDistMap = typename Network::template NodeMap<std::pair<DistType, NbHopsType>>; /**< Distance map. */
	using NodeDistMapSP = typename Network::template NodeMapSP<std::pair<DistType, NbHopsType>>; /**< Shared pointer to a distance map. */
	using EdgeLengthMap = typename NetSimlCalculator<Network, DistType, NbHopsType>::EdgeLengthMap; /**< Edge length map. */
	using EdgeLengthMapSP = typename NetSimlCalculator<Network, DistType, NbHopsType>::EdgeLengthMapSP; /**< Shared pointer to an edge length map. */

protected:
	using NetSimlCalculator<Network, DistType, NbHopsType>::selfSiml; /**< The value that should be assigned as similarity between a node and itself. */
	using NetSimlCalculator<Network, DistType, NbHopsType>::discNbHops; /**< The value that should be assigned as number of hops between disconnected nodes. */
	using NetSimlCalculator<Network, DistType, NbHopsType>::discDist; /**< The value that should be assigned as distance between disconnected nodes. */
	using NetSimlCalculator<Network, DistType, NbHopsType>::lambda; /**< Conductance. */
	Dijkstra<Network, DistType, NbHopsType> & dijkstra; /**< Dijkstra's algorithm. */
	std::shared_ptr<Network const> net; /**< The network. */
	EdgeLengthMapSP length; /**< The length map. */
	LengthMapIdType lengthMapId; /**< The ID of the length map. */
	CacheLevel cacheLevel; /**< The level of distance cache. */
	NodeDistMapSP nodeSimlCache; /**< Cache of distances from a node to all other nodes.*/
	std::queue<NodeID> cachedNodeIds; /**< Nodes in cache. */
	std::size_t maxNbNodesInCache; /**< Maximum number of nodes kept in cache (by default 10% of the nodes). */
	std::map<NodeID, NodeDistMapSP> netSimlCache; /**< Cache of similarity of a node to all other nodes.*/

private:
	double maxEdgeLength = 0; /**< Maximum length of any edge. */
	double avgEdgeLength = 0; /**< Average length of any edge. */

	/**
	 * Compute the similarity between two nodes given their similarity maps.
	 * @param i First node.
	 * @param j Second node.
	 * @param simi Similarity map of node i.
	 * @param simj Similarity map of node j.
	 * @return A tuple containing similarity, distance and number of hops between i and j.
	 */
	std::tuple<DistType, DistType, NbHopsType> getSiml(NodeID const & i,
			NodeID const & j, NodeDistMapSP const simi,
			NodeDistMapSP const simj) const;

	/**
	 * Compute the similarity between two nodes given their similarity maps.
	 * @param i First node.
	 * @param j Second node.
	 * @param simi Similarity map of node i.
	 * @param simj Similarity map of node j.
	 * @return A tuple containing similarity, distance and number of hops between i and j.
	 */
	std::tuple<DistType, DistType, NbHopsType> getDirSiml(NodeID const & i,
			NodeID const & j, NodeDistMapSP const simi) const;

public:
	/**
	 * Constructor.
	 * @param dijkstra A Dijkstra algorithm object.
	 * @param length The length map.
	 * @param cacheLevel The cache level.
	 */
	ESPSimlCalculator(Dijkstra<Network, DistType, NbHopsType> & dijkstra,
			EdgeLengthMapSP length, CacheLevel cacheLevel =
					CacheLevel::NetworkCache) :
			dijkstra(dijkstra), length(length), cacheLevel(cacheLevel) {
		net = dijkstra.getNet();
		lengthMapId = dijkstra.registerLengthMap(length);
		maxNbNodesInCache = net->getNbNodes() / 10;
		for (auto it = length->cbegin(); it != length->cend(); ++it) {
			if (it->second > maxEdgeLength) {
				maxEdgeLength = it->second;
			}
			avgEdgeLength += it->second;
		}
		avgEdgeLength /= length->size();
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	ESPSimlCalculator(ESPSimlCalculator const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	ESPSimlCalculator & operator =(ESPSimlCalculator const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	ESPSimlCalculator(ESPSimlCalculator && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	ESPSimlCalculator & operator =(ESPSimlCalculator && that) = default;

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The similarity between i and j.
	 */
	virtual std::tuple<DistType, DistType, NbHopsType> getSiml(
			NodeID const & i, NodeID const & j);

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The directed similarity between i and j.
	 */
	virtual std::tuple<DistType, DistType, NbHopsType> getDirSiml(
			NodeID const & i, NodeID const & j);

	/**
	 * @param e Edge.
	 * @return The similarity between the two ends of the edge.
	 */
	virtual std::tuple<DistType, DistType, NbHopsType> getSiml(
			Edge const & e) {
		return getSiml(net->start(e), net->end(e));
	}

	/**
	 * @param e An input edge.
	 * @return The directed similarity between start and end of e.
	 */
	virtual std::tuple<DistType, DistType, NbHopsType> getDirSiml(
			Edge const & e) {
		return getDirSiml(net->start(e), net->end(e));
	}

	/**
	 * @return The maximum number of nodes allowed in cache.
	 */
	std::size_t getMaxNbNodesInCache() const {
		return maxNbNodesInCache;
	}

	/**
	 * Set the maximum number of nodes in cache (for each node a map of similarity to all other nodes is kept in memory).
	 * @param maxNbNodesInCache New value for the maximum number of nodes in cache.
	 */
	void setMaxNbNodesInCache(std::size_t maxNbNodesInCache) {
		this->maxNbNodesInCache = maxNbNodesInCache;
	}

	/**
	 * Destructor.
	 */
	virtual ~ESPSimlCalculator() {
		dijkstra.unregisterLengthMap(lengthMapId);
	}

};

/**
 * @brief Exact shortest path distance calculator.
 * @details This class offers an additional layer over Dijkstra which provides
 * memory management functionalities by caching computed distances.
 * @tparam Network The network type.
 * @tparam DistType The distance type (can be an integer or floating point type).
 * @tparam NbHopsType The type of the number of hops (must be an integer type).
 */
template<typename Network = UNetwork<>, typename DistType = double,
		typename NbHopsType = std::size_t> class ESPLSimlCalculator: public NetSimlCalculator<
		Network, DistType, NbHopsType> {

public:
	using LengthMapIdType= typename NetSimlCalculator<Network, DistType, NbHopsType>::LengthMapIdType; /**< Length map ID type. */
	using NetworkSP = typename NetSimlCalculator<Network, DistType, NbHopsType>::NetworkSP; /**< Shared pointer to network. */
	using NodeID = typename NetSimlCalculator<Network, DistType, NbHopsType>::NodeID; /**< Nodes ID type. */
	using Label = typename NetSimlCalculator<Network, DistType, NbHopsType>::Label; /**< Nodes label type. */
	using Edge = typename NetSimlCalculator<Network, DistType, NbHopsType>::Edge; /**< Edge type. */
	using NodeSDistMap = typename Network::template NodeSMap<std::pair<DistType, NbHopsType>>; /**< Distance map. */
	using NodeSDistMapSP = typename Network::template NodeSMapSP<std::pair<DistType, NbHopsType>>; /**< Shared pointer to a distance map. */
	using EdgeLengthMap = typename NetSimlCalculator<Network, DistType, NbHopsType>::EdgeLengthMap; /**< Edge length map. */
	using EdgeLengthMapSP = typename NetSimlCalculator<Network, DistType, NbHopsType>::EdgeLengthMapSP; /**< Shared pointer to an edge length map. */

protected:
	using NetSimlCalculator<Network, DistType, NbHopsType>::selfSiml; /**< The value that should be assigned as similarity between a node and itself. */
	using NetSimlCalculator<Network, DistType, NbHopsType>::discNbHops; /**< The value that should be assigned as number of hops between disconnected nodes. */
	using NetSimlCalculator<Network, DistType, NbHopsType>::discDist; /**< The value that should be assigned as distance between disconnected nodes. */
	using NetSimlCalculator<Network, DistType, NbHopsType>::lambda; /**< Conductance. */
	Dijkstra<Network, DistType, NbHopsType> & dijkstra; /**< Dijkstra's algorithm. */
	EdgeLengthMapSP length; /**< The length map. */
	std::size_t lim; /**< Limit of number of hops. */
	std::shared_ptr<Network const> net; /**< The network. */
	LengthMapIdType lengthMapId; /**< The ID of the length map. */
	CacheLevel cacheLevel; /**< The level of distance cache. */
	NodeSDistMapSP nodeSimlCache; /**< Cache of distances from a node to all other nodes.*/
	std::queue<NodeID> cachedNodeIds; /**< Nodes in cache. */
	std::size_t maxNbNodesInCache; /**< Maximum number of nodes kept in cache (by default 10% of the nodes). */
	std::map<NodeID, NodeSDistMapSP> netSimlCache; /**< Cache of similarity of a node to all other nodes.*/
	bool useHops = true; /**< Use hops when computing similarity. If false, distance is used instead. */

private:

	double maxEdgeLength = 0; /**< Maximum length of any edge. */
	double avgEdgeLength = 0; /**< Average length of any edge. */

	/**
	 * Compute the similarity between two nodes given their similarity maps.
	 * @param i First node.
	 * @param j Second node.
	 * @param simi Similarity map of node i.
	 * @param simj Similarity map of node j.
	 * @return A tuple containing similarity, distance and number of hops between i and j.
	 */
	std::tuple<DistType, DistType, NbHopsType> getSiml(NodeID const & i,
			NodeID const & j, NodeSDistMapSP const simi,
			NodeSDistMapSP const simj) const;

	/**
	 * Compute the similarity between two nodes given their similarity maps.
	 * @param i First node.
	 * @param j Second node.
	 * @param simi Similarity map of node i.
	 * @param simj Similarity map of node j.
	 * @return A tuple containing similarity, distance and number of hops between i and j.
	 */
	std::tuple<DistType, DistType, NbHopsType> getDirSiml(NodeID const & i,
			NodeID const & j, NodeSDistMapSP const simi) const;

	/**
	 * @return The sparse distance from srcNode.
	 * @param i The source node.
	 */
	NodeSDistMapSP getSDistMap(NodeID const & i);

	inline double phi(double k) const {
		return 1 + std::log(k);
	}

public:
	/**
	 * Constructor.
	 * @param dijkstra A Dijkstra algorithm object.
	 * @param length The length map.
	 * @param lim Horizon limit.
	 * @param cacheLevel The cache level.
	 */
	ESPLSimlCalculator(Dijkstra<Network, DistType, NbHopsType> & dijkstra,
			EdgeLengthMapSP length, std::size_t lim, CacheLevel cacheLevel =
					CacheLevel::NetworkCache) :
			dijkstra(dijkstra), length(length), lim(lim), cacheLevel(cacheLevel) {
		net = dijkstra.getNet();
		lengthMapId = dijkstra.registerLengthMap(length);
		maxNbNodesInCache = net->getNbNodes() / 10;
		for (auto it = length->cbegin(); it != length->cend(); ++it) {
			if (it->second > maxEdgeLength) {
				maxEdgeLength = it->second;
			}
			avgEdgeLength += it->second;
		}
		avgEdgeLength /= length->size();
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	ESPLSimlCalculator(ESPLSimlCalculator const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	ESPLSimlCalculator & operator =(ESPLSimlCalculator const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	ESPLSimlCalculator(ESPLSimlCalculator && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	ESPLSimlCalculator & operator =(ESPLSimlCalculator && that) = default;

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The similarity between i and j.
	 */
	virtual std::tuple<DistType, DistType, NbHopsType> getSiml(
			NodeID const & i, NodeID const & j);

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The directed similarity between i and j.
	 */
	virtual std::tuple<DistType, DistType, NbHopsType> getDirSiml(
			NodeID const & i, NodeID const & j);

	/**
	 * @return A sparse similarity map of nodes having non-zero similarity to a given source node.
	 * @param srcNode The source node.
	 */
	virtual NodeSDistMapSP getNnzSimlMap(NodeID const & srcNode);

	/**
	 * @return A sparse similarity map of nodes having non-zero similarity to a given source node.
	 * Only nodes not connected to srcNode are considered. The node srcNode itself is also excluded.
	 * @param srcNode The source node.
	 */
	virtual NodeSDistMapSP getNnzSimlMapNoNeighb(NodeID const & srcNode);

	/**
	 * @param e An input edge.
	 * @return The similarity between start and end of e.
	 */
	virtual std::tuple<DistType, DistType, NbHopsType> getSiml(
			Edge const & e) {
		return getSiml(net->start(e), net->end(e));
	}

	/**
	 * @param e An input edge.
	 * @return The directed similarity between start and end of e.
	 */
	virtual std::tuple<DistType, DistType, NbHopsType> getDirSiml(
			Edge const & e) {
		return getDirSiml(net->start(e), net->end(e));
	}

	/**
	 * @return The maximum number of nodes allowed in cache.
	 */
	std::size_t getMaxNbNodesInCache() const {
		return maxNbNodesInCache;
	}

	/**
	 * Set the maximum number of nodes in cache (for each node a map of similarity to all other nodes is kept in memory).
	 * @param maxNbNodesInCache New value for the maximum number of nodes in cache.
	 */
	void setMaxNbNodesInCache(std::size_t maxNbNodesInCache) {
		this->maxNbNodesInCache = maxNbNodesInCache;
	}

	/**
	 * @return The limit on the number of hops.
	 */
	std::size_t getLim() const {
		return lim;
	}

	/**
	 * @return True if the number of hops is used to compute similarity.
	 */
	bool isUseHops() const {
		return useHops;
	}

	/**
	 * @param useHops Set whether the number of hops is used to compute similarity.
	 */
	void setUseHops(bool useHops) {
		this->useHops = useHops;
	}

	/**
	 * Destructor.
	 */
	virtual ~ ESPLSimlCalculator() {
		dijkstra.unregisterLengthMap(lengthMapId);
	}

};

/**
 * @brief Exact shortest path distance calculator.
 * @details This class offers an additional layer over Dijkstra which provides
 * memory management functionalities by caching computed distances.
 * @tparam Network The network type.
 * @tparam DistType The distance type (can be an integer or floating point type).
 * @tparam NbHopsType The type of the number of hops (must be an integer type).
 */
template<typename Network = UNetwork<>, typename DistType = double,
		typename NbHopsType = std::size_t> class ESPIndSimlCalculator: public NetIndSimlCalculator<
		Network, DistType, NbHopsType> {

public:
	using LengthMapIdType= typename NetIndSimlCalculator<Network, DistType, NbHopsType>::LengthMapIdType; /**< Length map ID type. */
	using NetworkSP = typename NetIndSimlCalculator<Network, DistType, NbHopsType>::NetworkSP; /**< Shared pointer to network. */
	using NodeID = typename NetIndSimlCalculator<Network, DistType, NbHopsType>::NodeID; /**< Nodes ID type. */
	using Label = typename NetIndSimlCalculator<Network, DistType, NbHopsType>::Label; /**< Nodes label type. */
	using Edge = typename NetIndSimlCalculator<Network, DistType, NbHopsType>::Edge; /**< Edge type. */
	using NodeDistMap = typename Network::template NodeMap<std::pair<DistType, NbHopsType>>; /**< Distance map. */
	using NodeDistMapSP = typename Network::template NodeMapSP<std::pair<DistType, NbHopsType>>; /**< Shared pointer to a distance map. */
	using EdgeLengthMap = typename NetIndSimlCalculator<Network, DistType, NbHopsType>::EdgeLengthMap; /**< Edge length map. */
	using EdgeLengthMapSP = typename NetIndSimlCalculator<Network, DistType, NbHopsType>::EdgeLengthMapSP; /**< Shared pointer to an edge length map. */

protected:
	using NetIndSimlCalculator<Network, DistType, NbHopsType>::selfSiml; /**< The value that should be assigned as similarity between a node and itself. */
	using NetIndSimlCalculator<Network, DistType, NbHopsType>::discNbHops; /**< The value that should be assigned as number of hops between disconnected nodes. */
	using NetIndSimlCalculator<Network, DistType, NbHopsType>::discDist; /**< The value that should be assigned as distance between disconnected nodes. */
	using NetIndSimlCalculator<Network, DistType, NbHopsType>::lambda; /**< Conductance. */
	Dijkstra<Network, DistType, NbHopsType> & dijkstra; /**< Dijkstra's algorithm. */
	std::shared_ptr<Network const> net; /**< The network. */
	EdgeLengthMapSP length; /**< The length map. */
	LengthMapIdType lengthMapId; /**< The ID of the length map. */
	std::map<std::pair<NodeID, NodeID>,
			std::vector<std::pair<DistType, NbHopsType>>> dirSimlCache; /**< Cache of directed similarity.*/

private:
	double maxEdgeLength = 0; /**< Maximum length of any edge. */
	double avgEdgeLength = 0; /**< Average length of any edge. */

	/**
	 * Compute the similarity between two nodes given their similarity maps.
	 * @param dikj Distance between i and j, followed by distances between i and all neighbors of j.
	 * @return A tuple containing similarity, distance and number of hops between i and j.
	 */
	std::tuple<DistType, DistType, NbHopsType> getDirIndSiml(
			std::vector<std::pair<DistType, NbHopsType>> const & dikj) const;

public:
	/**
	 * Constructor.
	 * @param dijkstra A Dijkstra algorithm object.
	 * @param length The length map.
	 */
	ESPIndSimlCalculator(Dijkstra<Network, DistType, NbHopsType> & dijkstra,
			EdgeLengthMapSP length) :
			dijkstra(dijkstra), length(length) {
		net = dijkstra.getNet();
		lengthMapId = dijkstra.registerLengthMap(length);
		for (auto it = length->cbegin(); it != length->cend(); ++it) {
			if (it->second > maxEdgeLength) {
				maxEdgeLength = it->second;
			}
			avgEdgeLength += it->second;
		}
		avgEdgeLength /= length->size();
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	ESPIndSimlCalculator(ESPIndSimlCalculator const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	ESPIndSimlCalculator & operator =(ESPIndSimlCalculator const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	ESPIndSimlCalculator(ESPIndSimlCalculator && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	ESPIndSimlCalculator & operator =(ESPIndSimlCalculator && that) = default;

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The similarity between i and j ignoring the edge between i and j.
	 */
	virtual std::tuple<DistType, DistType, NbHopsType> getIndSiml(
			NodeID const & i, NodeID const & j);

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The directed similarity between i and j ignoring the edge between i and j.
	 */
	virtual std::tuple<DistType, DistType, NbHopsType> getDirIndSiml(
			NodeID const & i, NodeID const & j);

	/**
	 * Destructor.
	 */
	virtual ~ESPIndSimlCalculator() {
		dijkstra.unregisterLengthMap(lengthMapId);
	}

};

/**
 * @brief Exact shortest path distance calculator on a directed network with limits on the number of hops.
 * @details This class offers an additional layer over Dijkstra which provides
 * memory management functionalities by caching computed distances.
 * @tparam Network The network type.
 * @tparam DistType The distance type (can be an integer or floating point type).
 * @tparam NbHopsType The type of the number of hops (must be an integer type).
 */
template<typename Network = DNetwork<>, typename DistType = double,
		typename NbHopsType = std::size_t> class DESPLDistCalculator: public NetDistCalculator<
		Network, DistType, NbHopsType> {

public:
	using LengthMapIdType= typename NetDistCalculator<Network, DistType, NbHopsType>::LengthMapIdType; /**< Length map ID type. */
	using NetworkSP = typename NetDistCalculator<Network, DistType, NbHopsType>::NetworkSP; /**< Shared pointer to network. */
	using NodeID = typename NetDistCalculator<Network, DistType, NbHopsType>::NodeID; /**< Nodes ID type. */
	using Label = typename NetDistCalculator<Network, DistType, NbHopsType>::Label; /**< Nodes label type. */
	using Edge = typename NetDistCalculator<Network, DistType, NbHopsType>::Edge; /**< Edge type. */
	using NodeDistMap = typename NetDistCalculator<Network, DistType, NbHopsType>::NodeDistMap; /**< Distance map. */
	using NodeDistMapSP = typename NetDistCalculator<Network, DistType, NbHopsType>::NodeDistMapSP; /**< Shared pointer to a distance map. */
	using NodeSDistMap = typename Network::template NodeSMap<std::pair<DistType, NbHopsType>>; /**< Distance map. */
	using NodeSDistMapSP = typename Network::template NodeSMapSP<std::pair<DistType, NbHopsType>>; /**< Shared pointer to a distance map. */
	using EdgeLengthMap = typename NetDistCalculator<Network, DistType, NbHopsType>::EdgeLengthMap; /**< Edge length map. */
	using EdgeLengthMapSP = typename NetDistCalculator<Network, DistType, NbHopsType>::EdgeLengthMapSP; /**< Shared pointer to an edge length map. */

protected:
	using NetDistCalculator<Network, DistType, NbHopsType>::discNbHops; /**< The value that should be assigned as number of hops between disconnected nodes. */
	using NetDistCalculator<Network, DistType, NbHopsType>::discDist; /**< The value that should be assigned as distance between disconnected nodes. */
	Dijkstra<Network, DistType, NbHopsType> & dijkstra; /**< Dijkstra's algorithm. */
	EdgeLengthMapSP length; /**< The length map. */
	std::size_t lim; /**< Limit of number of hops. */
	std::shared_ptr<Network const> net; /**< The network. */
	LengthMapIdType lengthMapId; /**< The ID of the length map. */
	CacheLevel cacheLevel; /**< The level of distance cache. */
	NodeSDistMapSP nodeDistCache; /**< Cache of distances from a node to all other nodes.*/
	std::queue<NodeID> cachedNodeIds; /**< Nodes in cache. */
	std::size_t maxNbNodesInCache; /**< Maximum number of nodes kept in cache (by default 10% of the nodes). */
	std::map<NodeID, NodeSDistMapSP> netDistCache; /**< Cache of distance of a node to all other nodes.*/

public:
	/**
	 * Constructor.
	 * @param dijkstra A Dijkstra algorithm object.
	 * @param lim Horizon limit.
	 * @param length The length map.
	 * @param cacheLevel The cache level.
	 */
	DESPLDistCalculator(Dijkstra<Network, DistType, NbHopsType> & dijkstra,
			EdgeLengthMapSP length, std::size_t lim, CacheLevel cacheLevel =
					CacheLevel::NetworkCache) :
			dijkstra(dijkstra), length(length), lim(lim), cacheLevel(cacheLevel) {
		net = dijkstra.getNet();
		lengthMapId = dijkstra.registerLengthMap(length);
		maxNbNodesInCache = net->getNbNodes() / 10;
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	DESPLDistCalculator(DESPLDistCalculator const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	DESPLDistCalculator & operator =(DESPLDistCalculator const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	DESPLDistCalculator(DESPLDistCalculator && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	DESPLDistCalculator & operator =(DESPLDistCalculator && that) = default;

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The distance between i and j.
	 */
	virtual std::pair<DistType, NbHopsType> getDist(NodeID const & i,
			NodeID const & j);

	/**
	 * @param i Index of the source node.
	 * @param j Index of the end node.
	 * @return The distance between i and j ignoring the edge between i and j.
	 */
	virtual std::pair<DistType, NbHopsType> getIndDist(NodeID const & i,
			NodeID const & j) {
		throw std::runtime_error("Not implemented");
	}

	/**
	 * @return The distance from srcNode.
	 * @param i The source node.
	 */
	virtual NodeDistMapSP getDist(NodeID const & i);

	/**
	 * @return The sparse distance from srcNode.
	 * @param i The source node.
	 */
	virtual NodeSDistMapSP getSDistMap(NodeID const & i);

	/**
	 * @return A sparse distance map of nodes having finite distance to a given source node.
	 * Only nodes not connected to srcNode are considered. The node srcNode itself is also excluded.
	 * @param srcNode The source node.
	 */
	virtual NodeSDistMapSP getFinDistMapNoNeighb(NodeID const & srcNode);

	/**
	 * @param e An input edge.
	 * @return The distance between start and end of e.
	 */
	virtual std::pair<DistType, NbHopsType> getDist(Edge const & e) {
		return getDist(net->start(e), net->end(e));
	}

	/**
	 * @return The maximum number of nodes allowed in cache.
	 */
	std::size_t getMaxNbNodesInCache() const {
		return maxNbNodesInCache;
	}

	/**
	 * Set the maximum number of nodes in cache (for each node a map of distance to all other nodes is kept in memory).
	 * @param maxNbNodesInCache New value for the maximum number of nodes in cache.
	 */
	void setMaxNbNodesInCache(std::size_t maxNbNodesInCache) {
		this->maxNbNodesInCache = maxNbNodesInCache;
	}

	/**
	 * @return The limit on the number of hops.
	 */
	std::size_t getLim() const {
		return lim;
	}

	/**
	 * Destructor.
	 */
	virtual ~ DESPLDistCalculator() {
		dijkstra.unregisterLengthMap(lengthMapId);
	}

};

} /* namespace LinkPred */

#endif /* DISTCALCULATOR_HPP_ */
