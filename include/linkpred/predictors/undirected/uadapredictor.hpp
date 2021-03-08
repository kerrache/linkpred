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
 * @brief Contains the implementation of an Adamic Adar index link predictor.
 */

#ifndef UADAPREDICTOR_HPP_
#define UADAPREDICTOR_HPP_

#include <linkpred/predictors/undirected/ulpredictor.hpp>
#include <memory>

namespace LinkPred {

/**
 * @brief Adamic Adar index link predictor.
 * @tparam Network The network type.
 * @tparam EdgeRndIt A random iterator type used to iterate on edges.
 * @tparam ScoreRndIt A random iterator type used to iterate on scores.
 */
template<typename Network = UNetwork<>,
		typename EdgeRndIt = typename std::vector<typename Network::Edge>::const_iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator,
		typename EdgeRndOutIt = typename std::vector<typename Network::Edge>::iterator> class UADAPredictor: public ULPredictor<
		Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt> {

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
	using NodeID = typename ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::NodeID; /**< The node IDs type. */
	using Edge = typename ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::Edge; /**< The edges type. */

public:
	/**
	 * @param net The network.
	 */
	UADAPredictor(std::shared_ptr<Network const> net) :
			ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>(net) {
		name = "ADA";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	UADAPredictor(UADAPredictor const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	UADAPredictor & operator =(UADAPredictor const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	UADAPredictor(UADAPredictor && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	UADAPredictor & operator =(UADAPredictor && that) = default;

	/**
	 * Initialize the predictor.
	 */
	virtual void init() {
	}

	/**
	 * Learning.
	 */
	virtual void learn() {
	}

	/**
	 * Compute the score of a single edge.
	 * @param e The edge.
	 * @return The score of e.
	 */
	virtual double score(Edge const & e) {
		auto srcNode = Network::start(e);
		auto endNode = Network::end(e);
		std::vector<NodeID> cn;
		net->getCommonNeighbors(srcNode, endNode, std::back_inserter(cn));

		double sum = 0;
		for (auto it = cn.begin(); it != cn.end(); ++it) {
			sum += 1.0 / (std::log(net->getDeg(*it)));
		}
		return sum;
	}

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
	 * Destructor.
	 */
	virtual ~UADAPredictor() = default;

};

} /* namespace LinkPred */

#endif /* UADAPREDICTOR_HPP_ */
