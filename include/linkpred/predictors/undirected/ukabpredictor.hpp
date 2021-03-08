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
 * @ingroup Predictors
 * @brief Contains the implementation of a scalable popularity-similarity link predictor.
 */

#ifndef UKABPREDICTOR_HPP_
#define UKABPREDICTOR_HPP_

#include "linkpred/predictors/undirected/ulpredictor.hpp"
#include "linkpred/graphalg/shortestpaths/netdistcalculator.hpp"
#include "linkpred/utils/log.hpp"
#include <memory>
#include <cmath>
#include <limits>
#include <iostream>

namespace LinkPred {

/**
 * @brief A scalable popularity similarity link predictor proposed in: "Kerrache, S., Alharbi, R. & Benhidour, H.
 * A Scalable Similarity-Popularity Link Prediction Method. Sci Rep 10, 6394 (2020)". https://doi.org/10.1038/s41598-020-62636-1.
 * @tparam Network The network type.
 * @tparam EdgeRndIt A random iterator type used to iterate on edges.
 * @tparam ScoreRndIt A random iterator type used to iterate on scores.
 */
template<typename Network = UNetwork<>,
		typename EdgeRndIt = typename std::vector<typename Network::Edge>::const_iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator,
		typename EdgeRndOutIt = typename std::vector<typename Network::Edge>::iterator> class UKABPredictor: public ULPredictor<
		Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt> {

public:
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::net; /**< The network. */
#ifdef LINKPRED_WITH_OPENMP
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt,
			EdgeRndOutIt>::parallel; /**< Whether the predictor runs in parallel. */
#endif
#ifdef LINKPRED_WITH_MPI
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt,
			EdgeRndOutIt>::comm; /**< The MPI communicator. */
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt,
			EdgeRndOutIt>::distributed; /**< Enable/disable distributed parallelism. */
#endif
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::name; /**< The name of the predictor. */
	using NodeID = typename ULPredictor<Network, EdgeRndIt,ScoreRndIt, EdgeRndOutIt>::NodeID; /**< The node IDs type. */
	using Edge = typename ULPredictor<Network, EdgeRndIt,ScoreRndIt, EdgeRndOutIt>::Edge; /**< The edges type. */

protected:

	CacheLevel cacheLevel = CacheLevel::NetworkCache; /**< The cache level. */
	std::size_t horizLim = 2; /**< Horizon limit. */

private:
	Dijkstra<Network, double, std::size_t> dijkstra; /**< Dijkstra's algorithm. */
	typename ESPLDistCalculator<Network, double, std::size_t>::EdgeLengthMapSP length; /**< The length map. */
	std::shared_ptr<ESPLDistCalculator<Network, double, std::size_t>> distCalc; /**< The distance calculator. */
	std::size_t maxDeg = 0; /**< The maximum degree in the network. */

	/**
	 * The function phi.
	 * @param k The input (typically, a node degree).
	 * @return Transformed input.
	 */
	inline double phi(double k) const {
		return std::log(k + 1);
	}

	/**
	 * Compute the normalized transformed sum of the degrees.
	 * @param ki Degree of node i.
	 * @param kj Degree of node j.
	 * @return The normalized transformed sum of the degrees.
	 */
	inline double sum(double ki, double kj) const {
		return (phi(ki) + phi(kj)) / (2 * phi(maxDeg));
	}

	/**
	 * Compute the neighborhood attraction force between two nodes.
	 * @param i Source node.
	 * @param j End node.
	 * @param ki Degree of i.
	 * @param kj Degree of j.
	 * @return the neighborhood attraction force between i and j.
	 */
	inline double nsc(NodeID i, NodeID j, double ki, double kj) const {
		std::vector<NodeID> cn;
		net->getCommonNeighbors(i, j, std::back_inserter(cn));
		double ns = 1;
		for (auto it = cn.begin(); it != cn.end(); ++it) {
			double kk = net->getDeg(*it);
			ns *= phi(kk) / phi(maxDeg);
		}
		return 1 - ns;
	}

	/**
	 * Compute the length of the edge.
	 * @param edge The edge.
	 * @return The edge length using PAT method.
	 */
	inline double edgeLength(Edge const &edge) const {
		auto i = Network::start(edge);
		auto j = Network::end(edge);

		double ki = net->getDeg(i);
		double kj = net->getDeg(j);
		double etaij = nsc(i, j, ki, kj);
		double popij = sum(ki, kj);
		double sij = 2 * popij / (1 + etaij);
		return sij;
	}

	/**
	 * Create the length map from network topology.
	 */
	void createLengthMap();

	/**
	 * @return A sparse node map containing the distance to all nodes having a nonzero similarity with i.
	 * @param i Source node.
	 */
	typename Network::template NodeSMapSP<std::pair<double, size_t>> getFiniteDistMap(
			NodeID i) {
		auto distMap = distCalc->getFinDistMapNoNeighb(i);
		for (auto it = distMap->begin(); it != distMap->end(); ++it) {
			it->second.first = 1.0 + it->second.first;
		}
		return distMap;
	}

	/**
	 * Computes the similarity between two nodes.
	 * @param i Source node.
	 * @param j End node.
	 */
	double getDist(NodeID i, NodeID j) {
		double dist = 1 + distCalc->getDist(i, j).first;
		return dist;
	}

	/**
	 * Computes the similarity between the two ends of an edge.
	 * @param e The edge.
	 * @return The similarity between the two ends of e.
	 */
	double getDist(Edge const &e) {
		return getDist(net->start(e), net->end(e));
	}

	/**
	 * Compute the score of a single edge.
	 * @param e The edge.
	 * @param dist Distance associated with e.
	 * @return The score of e.
	 */
	double score(Edge const &e, double dist);

public:
	/**
	 * @param net The network.
	 */
	UKABPredictor(std::shared_ptr<Network const> net) :
			ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>(net), dijkstra(
					net) {
		name = "KAB";
		maxDeg = net->getMaxDeg();
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	UKABPredictor(UKABPredictor const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	UKABPredictor& operator =(UKABPredictor const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	UKABPredictor(UKABPredictor &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	UKABPredictor& operator =(UKABPredictor &&that) = default;

	/**
	 * Initialize the predictor.
	 */
	virtual void init();

	/**
	 * Learn.
	 */
	virtual void learn();

	/**
	 * Compute the score of a single edge.
	 * @param e The edge.
	 * @return The score of e.
	 */
	virtual double score(Edge const &e);

	/**
	 * Predict the links.
	 * @param begin Beginning of the links to be predicted.
	 * @param end end of the links to be predicted.
	 * @param scores Beginning of scores.
	 */
	virtual void predict(EdgeRndIt begin, EdgeRndIt end, ScoreRndIt scores);

	/**
	 * Finds the k negative edges with the top score. Ties are broken randomly.
	 * @param k The number of edges to find.
	 * @param eit An output iterator where the edges are written.
	 * @param sit An output iterator where the scores are written. The scores are written in the same order as the edges.
	 * @return The number of negative edges inserted. It is the minimum between k and the number of negative edges in the network.
	 */
	virtual std::size_t top(std::size_t k, EdgeRndOutIt eit, ScoreRndIt sit);

	/**
	 * @return The distances cache level.
	 */
	CacheLevel getCacheLevel() const {
		return cacheLevel;
	}

	/**
	 * Set the distances cache level.
	 * @param cacheLevel  The new distances cache level.
	 */
	void setCacheLevel(CacheLevel cacheLevel) {
		this->cacheLevel = cacheLevel;
	}

	/**
	 * @return The horizon limit.
	 */
	std::size_t getHorizLim() const {
		return horizLim;
	}

	/**
	 * @param horizLim New horizon limit.
	 */
	void setHorizLim(std::size_t horizLim) {
		this->horizLim = horizLim;
	}

	/**
	 * Destructor.
	 */
	virtual ~UKABPredictor() = default;

};

} /* namespace LinkPred */

#endif /* UKABPREDICTOR_HPP_ */
