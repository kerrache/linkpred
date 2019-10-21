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
 * @brief Contains the implementation of a link prediction algorithm that combines node sum (popularity) and local attraction.
 */

#ifndef USNSPREDICTOR_HPP_
#define USNSPREDICTOR_HPP_

#include <linkpred/predictors/ulpredictor.hpp>
#include "linkpred/utils/log.hpp"
#include <memory>
#include <cmath>
#include <limits>

namespace LinkPred {

/**
 * @brief A neighbors degree link predictor.
 * @tparam NetworkT The network type.
 * @tparam EdgesRandomIteratorT A random iterator type used to iterate on edges.
 * @tparam ScoresRandomIteratorT A random iterator type used to iterate on scores.
 */
template<typename NetworkT = UNetwork<>,
		typename EdgesRandomIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::const_iterator,
		typename ScoresRandomIteratorT = typename std::vector<double>::iterator,
		typename EdgesRandomOutputIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::iterator> class USNSPredictor: public ULPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT> {

public:
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::net;
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
			EdgesRandomOutputIteratorT>::name;
	using NodeIdType = typename ULPredictor<NetworkT, EdgesRandomIteratorT,ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::NodeIdType;
	using EdgeType = typename ULPredictor<NetworkT, EdgesRandomIteratorT,ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::EdgeType;

protected:

	std::size_t maxDeg = 0; /**< The maximum degree in the network. */
	double lambda = 0.5;

	inline double phi(double k) const {
		return std::log(k + 1);
	}

	inline double sum(double ki, double kj) const {
		return (phi(ki) + phi(kj)) / (2 * phi(maxDeg));
	}

	inline double nsc(NodeIdType srcNode, NodeIdType endNode, double ki,
			double kj) const {
		std::vector<NodeIdType> cn;
		net->getCommonNeighbors(srcNode, endNode, std::back_inserter(cn));
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
	 * @param seed The random number generator's seed.
	 */
	USNSPredictor(std::shared_ptr<NetworkT const> net) :
			ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
					EdgesRandomOutputIteratorT>(net) {
		name = "SNS";
		maxDeg = net->getMaxDeg();
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	USNSPredictor(USNSPredictor const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	USNSPredictor & operator =(USNSPredictor const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	USNSPredictor(USNSPredictor && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	USNSPredictor & operator =(USNSPredictor && that) = default;

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
	virtual double score(EdgeType const & e);

	/**
	 * Predict the links.
	 * @param begin Beginning of the links to be predicted.
	 * @param end end of the links to be predicted.
	 * @param scores Beginning of scores.
	 */
	virtual void predict(EdgesRandomIteratorT begin, EdgesRandomIteratorT end,
			ScoresRandomIteratorT scores);

	double getLambda() const {
		return lambda;
	}

	void setLambda(double lambda) {
		this->lambda = lambda;
	}

	/**
	 * Destructor.
	 */
	virtual ~USNSPredictor() = default;

};
}
/* namespace LinkPred */

#endif /* USNSPREDICTOR_HPP_ */
