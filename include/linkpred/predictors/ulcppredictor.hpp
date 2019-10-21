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
 * @brief Contains the implementation of a local path link predictor.
 */

#ifndef ULCPPREDICTOR_HPP_
#define ULCPPREDICTOR_HPP_

#include <linkpred/predictors/ulpredictor.hpp>
#include <memory>

namespace LinkPred {

/**
 * @brief Local path link predictor.
 * @tparam NetworkT The network type.
 * @tparam EdgesRandomIteratorT A random iterator type used to iterate on edges.
 * @tparam ScoresRandomIteratorT A random iterator type used to iterate on scores.
 */
template<typename NetworkT = UNetwork<>,
		typename EdgesRandomIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::const_iterator,
		typename ScoresRandomIteratorT = typename std::vector<double>::iterator,
		typename EdgesRandomOutputIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::iterator> class ULCPPredictor: public ULPredictor<
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

protected:
	double epsilon = 0.001; /**< The weight of paths of length 3. */

public:
	/**
	 * @param net The network.
	 */
	ULCPPredictor(std::shared_ptr<NetworkT const> net) :
			ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
					EdgesRandomOutputIteratorT>(net) {
		name = "LCP";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	ULCPPredictor(ULCPPredictor const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	ULCPPredictor & operator =(ULCPPredictor const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	ULCPPredictor(ULCPPredictor && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	ULCPPredictor & operator =(ULCPPredictor && that) = default;

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
	 * @return epsilon, the weight of paths of length 3.
	 */
	double getEpsilon() const {
		return epsilon;
	}

	/**
	 * Set epsilon, the weight of paths of length 3.
	 * @param epsilon The new weight of paths of length 3.
	 */
	void setEpsilon(double epsilon) {
		this->epsilon = epsilon;
	}

	/**
	 * Destructor.
	 */
	virtual ~ULCPPredictor() = default;
};
}
/* namespace LinkPred */

#endif /* ULCPPREDICTOR_HPP_ */
