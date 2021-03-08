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
 * @brief Contains the implementation of a shortest path link predictor.
 */

#ifndef USHPPREDICTOR_HPP_
#define USHPPREDICTOR_HPP_

#include <linkpred/predictors/undirected/ulpredictor.hpp>
#include "linkpred/graphalg/shortestpaths/netdistcalculator.hpp"
#include "linkpred/utils/log.hpp"
#include <memory>
#include <cmath>

namespace LinkPred {

/**
 * @brief A shortest path link predictor link predictor.
 * @tparam Network The network type.
 * @tparam EdgeRndIt A random iterator type used to iterate on edges.
 * @tparam ScoreRndIt A random iterator type used to iterate on scores.
 */
template<typename Network = UNetwork<>,
		typename EdgeRndIt = typename std::vector<typename Network::Edge>::const_iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator,
		typename EdgeRndOutIt = typename std::vector<typename Network::Edge>::iterator> class USHPPredictor: public ULPredictor<
		Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt> {

	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::net; /**< The network. */
#ifdef LINKPRED_WITH_OPENMP
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::parallel; /**< Whether the predictor runs in parallel. */
#endif
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::name; /**< The name of the predictor. */
	using NodeID = typename ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::NodeID; /**< The node IDs type. */
	using Edge = typename ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::Edge; /**< The edges type. */

public:
	/**
	 * @brief An enumeration of the different landmark positioning strategies.
	 */
	enum LandmarkStrategy {
		Random, /**< Landmarks are chosen randomly. */
		Hub, /**< The nodes with the highest degree are chosen. */
		IHub /**< The nodes with the lowest degree are chosen. */
	};

protected:
	long int seed = 0; /**< The random number generator seed. */
	RandomGen rng; /**< The random number generator. */
	Dijkstra<Network, double, std::size_t> dijkstra; /**< Dijkstra's algorithm. */
	typename NetDistCalculator<Network, double, std::size_t>::EdgeLengthMapSP length; /**< The length map. */
	std::shared_ptr<NetDistCalculator<Network, double, std::size_t>> distCalc; /**< The distance calculator. */
	CacheLevel cacheLevel = CacheLevel::NodeCache; /**< Cache level for network distances. */
	bool asp = false; /**< Whether to use approximate shortest paths instead of exact ones. */
	double landmarkRatio = 0; /**< The ratio of nodes used as landmarks if asp is on. */
	LandmarkStrategy landmarkStrategy = Random; /**< Landmark positioning strategy if asp is on. */
	std::set<NodeID> landmarks; /**< The landmarks. */

	/**
	 * Create the length map from network topology.
	 */
	void createLengthMap();

	/**
	 * Position the landmarks.
	 */
	void setLandmarks();

public:
	/**
	 * @param net The network.
	 * @param seed The random number generator's seed.
	 */
	USHPPredictor(std::shared_ptr<Network const> net, long int seed) :
			ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>(net), seed(
					seed), rng(seed), dijkstra(net) {
		name = "SHP";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	USHPPredictor(USHPPredictor const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	USHPPredictor & operator =(USHPPredictor const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	USHPPredictor(USHPPredictor && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	USHPPredictor & operator =(USHPPredictor && that) = default;

	/**
	 * Initialize the solver.
	 */
	virtual void init();

	/**
	 * Learn.
	 */
	virtual void learn();

	/**
	 * Predict the links.
	 * @param begin Beginning of the links to be predicted.
	 * @param end end of the links to be predicted.
	 * @param scores Beginning of scores.
	 */
	virtual void predict(EdgeRndIt begin, EdgeRndIt end, ScoreRndIt scores);

	/**
	 * Compute the score of a single edge.
	 * @param e The edge.
	 * @return The score of e.
	 */
	virtual double score(Edge const & e);

	/**
	 * @return Whether approximate shortest path distances are used.
	 */
	bool getAsp() const {
		return asp;
	}

	/**
	 * Set Whether approximate shortest path distances are used.
	 * @param asp Whether approximate shortest path distances should be used.
	 */
	void setAsp(bool asp) {
		this->asp = asp;
	}

	/**
	 * @return The landmark ratio.
	 */
	double getLandmarkRatio() const {
		return landmarkRatio;
	}

	/**
	 * Set the landmark ratio.
	 * @param landmarkRatio The new landmark ratio.
	 */
	void setLandmarkRatio(double landmarkRatio) {
		if (landmarkRatio < 0 || landmarkRatio > 1) {
			throw std::invalid_argument(
					"Invalid landmark ratio, must be between 0 and 1 inclusive");
		}
		this->landmarkRatio = landmarkRatio;
	}

	/**
	 * @return The landmark strategy.
	 */
	LandmarkStrategy getLandmarkStrategy() const {
		return landmarkStrategy;
	}

	/**
	 * Set the landmark strategy.
	 * @param landmarkStrategy The new landmark strategy.
	 */
	void setLandmarkStrategy(LandmarkStrategy landmarkStrategy) {
		this->landmarkStrategy = landmarkStrategy;
	}

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
	 * Destructor.
	 */
	virtual ~USHPPredictor() = default;

};
}
/* namespace LinkPred */

#endif /* USHPPREDICTOR_HPP_ */
