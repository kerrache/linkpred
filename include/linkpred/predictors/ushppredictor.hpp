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
 * @brief Contains the implementation of a shortest path link predictor.
 */

#ifndef USHPPREDICTOR_HPP_
#define USHPPREDICTOR_HPP_

#include <linkpred/predictors/ulpredictor.hpp>
#include "linkpred/core/netdistcalculator.hpp"
#include "linkpred/utils/log.hpp"
#include <memory>
#include <cmath>

namespace LinkPred {
/**
 * @brief A shortest path link predictor link predictor.
 * @tparam NetworkT The network type.
 * @tparam EdgesRandomIteratorT A random iterator type used to iterate on edges.
 * @tparam ScoresRandomIteratorT A random iterator type used to iterate on scores.
 */
template<typename NetworkT = UNetwork<>,
		typename EdgesRandomIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::const_iterator,
		typename ScoresRandomIteratorT = typename std::vector<double>::iterator,
		typename EdgesRandomOutputIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::iterator> class USHPPredictor: public ULPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT> {

	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::net; /**< The network. */
#ifdef WITH_OPENMP
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::parallel; /**< Whether the predictor runs in parallel. */
#endif
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::name; /**< The name of the predictor. */
	using NodeIdType = typename ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::NodeIdType; /**< The node IDs type. */
	using EdgeType = typename ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::EdgeType; /**< The edges type. */

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
	Dijkstra<NetworkT, double, std::size_t> dijkstra; /**< Dijkstra's algorithm. */
	typename NetDistCalculator<NetworkT, double, std::size_t>::EdgeLengthMapSP length; /**< The length map. */
	std::shared_ptr<NetDistCalculator<NetworkT, double, std::size_t>> distCalc; /**< The distance calculator. */
	CacheLevel cacheLevel = NodeCache; /**< Cache level for network distances. */
	bool asp = false; /**< Whether to use approximate shortest paths instead of exact ones. */
	double landmarkRatio = 0; /**< The ratio of nodes used as landmarks if asp is on. */
	LandmarkStrategy landmarkStrategy = Random; /**< Landmark positioning strategy if asp is on. */
	std::set<NodeIdType> landmarks; /**< The landmarks. */

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
	USHPPredictor(std::shared_ptr<NetworkT const> net, long int seed) :
			ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
					EdgesRandomOutputIteratorT>(net), seed(seed), rng(seed), dijkstra(
					net) {
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
	virtual void predict(EdgesRandomIteratorT begin, EdgesRandomIteratorT end,
			ScoresRandomIteratorT scores);

	/**
	 * Compute the score of a single edge.
	 * @param e The edge.
	 * @return The score of e.
	 */
	virtual double score(EdgeType const & e);

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
