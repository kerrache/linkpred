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
 * @brief Contains the implementation of a random link predictor.
 */

#ifndef URNDPREDICTOR_HPP_
#define URNDPREDICTOR_HPP_

#include <linkpred/predictors/ulpredictor.hpp>
#include "linkpred/utils/randomgen.hpp"
#include <memory>

namespace LinkPred {

/**
 * @brief A random link predictor.
 * @tparam NetworkT The network type.
 * @tparam EdgesRandomIteratorT A random iterator type used to iterate on edges.
 * @tparam ScoresRandomIteratorT A random iterator type used to iterate on scores.
 */
template<typename NetworkT = UNetwork<>,
		typename EdgesRandomIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::const_iterator,
		typename ScoresRandomIteratorT = typename std::vector<double>::iterator,
		typename EdgesRandomOutputIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::iterator> class URNDPredictor: public ULPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT> {

	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::net; /**< The network. */
#ifdef WITH_OPENMP
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::parallel; /**< Whether the predictor runs in parallel. */
#endif
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::name; /**< The name of the predictor. */
	using NodeIdType = typename ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::NodeIdType; /**< The node IDs type. */
	using EdgeType = typename ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::EdgeType; /**< The edges type. */

protected:
	RandomGen rng; /**< The random number generator. */

public:
	/**
	 * @param net The network.
	 * @param seed The random number generator's seed.
	 */
	URNDPredictor(std::shared_ptr<NetworkT const> net, long int seed) :
			ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
					EdgesRandomOutputIteratorT>(net), rng(seed) {
		name = "RND";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	URNDPredictor(URNDPredictor const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	URNDPredictor & operator =(URNDPredictor const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	URNDPredictor(URNDPredictor && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	URNDPredictor & operator =(URNDPredictor && that) = default;

	/**
	 * Initialize the solver.
	 */
	virtual void init() {
	}

	/**
	 * Learn.
	 */
	virtual void learn() {
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
	 * Compute the score of a single edge.
	 * @param e The edge.
	 * @return The score of e.
	 */
	virtual double score(EdgeType const & e);

	/**
	 * Destructor.
	 */
	virtual ~URNDPredictor() = default;
};

}
/* namespace LinkPred */

#endif /* URNDPREDICTOR_HPP_ */
