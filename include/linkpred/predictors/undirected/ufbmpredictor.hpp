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
 * @brief Contains the implementation of the fast blocking model link predictor.
 */

#ifndef UFBMPREDICTOR_HPP_
#define UFBMPREDICTOR_HPP_

#include <linkpred/predictors/undirected/ulpredictor.hpp>
#include "linkpred/numerical/linear/gfmatrix.hpp"
#include "linkpred/utils/randomgen.hpp"

namespace LinkPred {

/**
 * @brief Fast blocking model link predictor.
 * @details This is a C++ translation of the Matlab code provided by the authors.
 * @tparam Network The network type.
 * @tparam EdgeRndIt A random iterator type used to iterate on edges.
 * @tparam ScoreRndIt A random iterator type used to iterate on scores.
 */
template<typename Network = UNetwork<>,
		typename EdgeRndIt = typename std::vector<typename Network::Edge>::const_iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator,
		typename EdgeRndOutIt = typename std::vector<typename Network::Edge>::iterator> class UFBMPredictor: public ULPredictor<
		Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt> {

	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::net; /**< The network. */
#ifdef LINKPRED_WITH_OPENMP
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::parallel; /**< Whether the predictor runs in parallel. */
#endif
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::name; /**< The name of the predictor. */
	using NodeID = typename ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::NodeID; /**< The node IDs type. */
	using Edge = typename ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::Edge; /**< The edges type. */

protected:
	std::shared_ptr<GFMatrix> adj; /**< Adjacency matrix. */
	std::map<Edge, double> scores; /**< Links scores. */
	std::size_t maxIter = 50; /**< Number of samplings. */
	long int seed; /**< Random number generator seed. */
	RandomGen rng; /**< Random number generator. */

	/**
	 * @param train Adjacency matrix of the observed network.
	 * @param i The ith node in the network.
	 * @return A clique that includes node i.
	 */
	std::vector<std::size_t> cliqueFind(GFMatrix const & train,
			std::size_t i);

	/**
	 * @param s1: the node vector of the current network.
	 * @param train: adjacency matrix of the current network.
	 * @return A pair containing: (1): an obtained node vector which denotes a clique. (2): a node vector which denotes the remainder of the network.
	 */
	std::pair<std::vector<std::size_t>, std::vector<std::size_t>> blockG(
			std::vector<std::size_t> const & s1,
			GFMatrix const & train);

	/**
	 * @param train: adjacency matrix of the observed network.
	 * @return A pair containing: (1): The obtained communities. (2): The index of the first special community.
	 */
	std::pair<GFMatrix, std::size_t> divisionG(
			GFMatrix const & train);

	/**
	 * Main FBM loop.
	 * @param adj The adjacency matrix.
	 */
	void fbm(std::shared_ptr<GFMatrix> adj);

public:

	/**
	 * @param net The network.
	 * @param seed Random number generator seed.
	 */
	UFBMPredictor(std::shared_ptr<Network const> net, long int seed) :
			ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>(net), seed(
					seed), rng(seed) {
		name = "FBM";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	UFBMPredictor(UFBMPredictor const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	UFBMPredictor & operator =(UFBMPredictor const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	UFBMPredictor(UFBMPredictor && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	UFBMPredictor & operator =(UFBMPredictor && that) = default;

	/**
	 * Initialize the solver.
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
	virtual double score(Edge const & e) {
		return scores.at(e);
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
	virtual ~UFBMPredictor() = default;
};

}
/* namespace LinkPred */

#endif /* UFBMPREDICTOR_HPP_ */
