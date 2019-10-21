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
 * @brief Contains the implementation of an Adamic Adar index link predictor.
 */

#ifndef UADAPREDICTOR_HPP_
#define UADAPREDICTOR_HPP_

#include <linkpred/predictors/ulpredictor.hpp>
#include <memory>

namespace LinkPred {

/**
 * @brief Adamic Adar index link predictor.
 * @tparam NetworkT The network type.
 * @tparam EdgesRandomIteratorT A random iterator type used to iterate on edges.
 * @tparam ScoresRandomIteratorT A random iterator type used to iterate on scores.
 */
template<typename NetworkT = UNetwork<>,
		typename EdgesRandomIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::const_iterator,
		typename ScoresRandomIteratorT = typename std::vector<double>::iterator,
		typename EdgesRandomOutputIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::iterator> class UADAPredictor: public ULPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT> {

	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::net; /**< The network. */
#ifdef WITH_OPENMP
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::parallel; /**< Whether the predictor runs in parallel. */
#endif
#ifdef WITH_MPI
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::comm; /**< The MPI communicator. */
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::distributed; /**< Enable/disable distributed parallelism. */

#endif
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::name; /**< The name of the predictor. */
	using NodeIdType = typename ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::NodeIdType; /**< The node IDs type. */
	using EdgeType = typename ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::EdgeType; /**< The edges type. */

public:
	/**
	 * @param net The network.
	 */
	UADAPredictor(std::shared_ptr<NetworkT const> net) :
			ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
					EdgesRandomOutputIteratorT>(net) {
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
	virtual double score(EdgeType const & e) {
		auto srcNode = NetworkT::start(e);
		auto endNode = NetworkT::end(e);
		std::vector<NodeIdType> cn;
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
	virtual void predict(EdgesRandomIteratorT begin, EdgesRandomIteratorT end,
			ScoresRandomIteratorT scores);

	/**
	 * Finds the k negative edges with the top score. Ties are broken randomly.
	 * @param k The number of edges to find.
	 * @param eit An output iterator where the edges are written.
	 * @param sit An output iterator where the scores are written. The scores are written in the same order as the edges.
	 * @return The number of negative edges inserted. It is the minimum between k and the number of negative edges in the network.
	 */
	virtual std::size_t top(std::size_t k, EdgesRandomOutputIteratorT eit,
			ScoresRandomIteratorT sit);

	/**
	 * Destructor.
	 */
	virtual ~UADAPredictor() = default;

};
}
/* namespace LinkPred */

#endif /* UADAPREDICTOR_HPP_ */
