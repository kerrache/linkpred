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
 * @brief Contains the implementation of a link prediction algorithm based on the degrees of neighbors.
 */

#ifndef UNEDPREDICTOR_HPP_
#define UNEDPREDICTOR_HPP_

#include <linkpred/predictors/undirected/ulpredictor.hpp>
#include "linkpred/utils/log.hpp"
#include <memory>
#include <cmath>
#include <limits>

namespace LinkPred {

/**
 * @brief A neighbors degree link predictor.
 * @tparam Network The network type.
 * @tparam EdgeRndIt A random iterator type used to iterate on edges.
 * @tparam ScoreRndIt A random iterator type used to iterate on scores.
 */
template<typename Network = UNetwork<>,
		typename EdgeRndIt = typename std::vector<typename Network::Edge>::const_iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator,
		typename EdgeRndOutIt = typename std::vector<typename Network::Edge>::iterator> class UNEDPredictor: public ULPredictor<
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

private:

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

public:

	/**
	 * @param net The network.
	 */
	UNEDPredictor(std::shared_ptr<Network const> net) :
			ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>(net) {
		name = "NED";
		maxDeg = net->getMaxDeg();
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	UNEDPredictor(UNEDPredictor const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	UNEDPredictor& operator =(UNEDPredictor const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	UNEDPredictor(UNEDPredictor &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	UNEDPredictor& operator =(UNEDPredictor &&that) = default;

	/**
	 * Initialize the predictor.
	 */
	virtual void init() {
	}

	/**
	 * Learn.
	 */
	virtual void learn() {
	}

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
	 * Finds the k negative edges with the top score. Ties are broken randmly.
	 * @param k The number of edges to find.
	 * @param eit An output iterator where the edges are written.
	 * @param sit An output iterator where the scores are written. The scores are written in the same order as the edges.
	 * @return The number of negative edges inserted. It is the minimum between k and the number of negative edges in the network.
	 */
	virtual std::size_t top(std::size_t k, EdgeRndOutIt eit, ScoreRndIt sit);

	/**
	 * Destructor.
	 */
	virtual ~UNEDPredictor() = default;
};

} /* namespace LinkPred */

#endif /* UNEDPREDICTOR_HPP_ */
