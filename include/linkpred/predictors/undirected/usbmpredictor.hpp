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
 * @brief Contains the implementation of the stochastic block model link predictor.
 */

#ifndef USBMPREDICTOR_HPP_
#define USBMPREDICTOR_HPP_

#include <linkpred/predictors/undirected/ulpredictor.hpp>
#include <linkpred/predictors/undirected/usbmpredictor/sbm.hpp>

namespace LinkPred {

/**
 * @brief The stochastic block model link predictor.
 * @tparam Network The network type.
 * @tparam EdgeRndIt A random iterator type used to iterate on edges.
 * @tparam ScoreRndIt A random iterator type used to iterate on scores.
 */
template<typename Network = UNetwork<>,
		typename EdgeRndIt = typename std::vector<typename Network::Edge>::const_iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator,
		typename EdgeRndOutIt = typename std::vector<typename Network::Edge>::iterator> class USBMPredictor: public ULPredictor<
		Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt> {

	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::net; /**< The network. */
#ifdef LINKPRED_WITH_OPENMP
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::parallel; /**< Whether the predictor runs in parallel. */
#endif
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::name; /**< The name of the predictor. */
	using NodeID = typename ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::NodeID; /**< The node IDs type. */
	using Edge = typename ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::Edge; /**< The edges type. */

protected:
	long int seed = 0; /**< The random number generator seed. */
	RandomGen rng; /**< The random number generator. */
	struct node_gra * rGNet = nullptr; /**< Network in RGraph format. */
	std::size_t maxIter = 1000; /**< Maximum number of iterations. */
	double** scores = nullptr; /**< Links scores. */
	std::vector<struct node_gra *> nodesMap; /**< Map of node IDs. */

	/**
	 * Convert the network to an SBM graph structure.
	 * @param net The network.
	 * @return an SBM graph.
	 */
	struct node_gra * toRGraphNet(std::shared_ptr<Network const> net);

public:

	/**
	 * @param net The network.
	 * @param seed The random number generator's seed.
	 */
	USBMPredictor(std::shared_ptr<Network const> net, long int seed) :
			ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>(net), seed(
					seed), rng(seed) {
		name = "SBM";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	USBMPredictor(USBMPredictor const & that) = delete;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	USBMPredictor & operator =(USBMPredictor const & that) = delete;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	USBMPredictor(USBMPredictor && that) = delete;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	USBMPredictor & operator =(USBMPredictor && that) = delete;

	/**
	 * Initialize the solver.
	 */
	virtual void init();

	/**
	 * Learning.
	 */
	virtual void learn();

	/**
	 * Compute the score of a single edge.
	 * @param e The edge.
	 * @return The score of e.
	 */
	virtual double score(Edge const & e) {
		auto i = Network::start(e);
		auto j = Network::end(e);
		return scores[nodesMap[i]->num][nodesMap[j]->num];
	}

	/**
	 * Predict the links.
	 * @param begin Beginning of the links to be predicted.
	 * @param end end of the links to be predicted.
	 * @param scores Beginning of scores.
	 */
	virtual void predict(EdgeRndIt begin, EdgeRndIt end, ScoreRndIt scores);

	/**
	 * @return The maximum number of iterations.
	 */
	std::size_t getMaxIter() const {
		return maxIter;
	}

	/**
	 * Set the maximum number of iterations.
	 * @param maxIter The new maximum number of iterations.
	 */
	void setMaxIter(std::size_t maxIter) {
		this->maxIter = maxIter;
	}

	/**
	 * Destructor.
	 */
	virtual ~USBMPredictor();
};

} /* namespace LinkPred */

#endif /* USBMPREDICTOR_HPP_ */
