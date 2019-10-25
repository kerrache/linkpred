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
 * @brief Contains the implementation of the stochastic block model link predictor.
 */

#ifndef USBMPREDICTOR_HPP_
#define USBMPREDICTOR_HPP_

#include <linkpred/predictors/ulpredictor.hpp>
#include <linkpred/predictors/usbmpredictor/sbm.hpp>

namespace LinkPred {

/**
 * @brief The stochastic block model link predictor.
 * @tparam NetworkT The network type.
 * @tparam EdgesRandomIteratorT A random iterator type used to iterate on edges.
 * @tparam ScoresRandomIteratorT A random iterator type used to iterate on scores.
 */
template<typename NetworkT = UNetwork<>,
		typename EdgesRandomIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::const_iterator,
		typename ScoresRandomIteratorT = typename std::vector<double>::iterator,
		typename EdgesRandomOutputIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::iterator> class USBMPredictor: public ULPredictor<
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
	long int seed = 0; /**< The random number generator seed. */
	RandomGen rng; /**< The random number generator. */
	struct node_gra * rGNet = nullptr; /**< Network in RGraph format. */
	std::size_t maxIter = 10000; /**< Maximum number of iterations. */
	double** scores = nullptr; /**< Links scores. */
	std::vector<struct node_gra *> nodesMap; /**< Map of node IDs. */

	/**
	 * Convert the network to an SBM graph structure.
	 * @param net The network.
	 * @return an SBM graph.
	 */
	struct node_gra * toRGraphNet(std::shared_ptr<NetworkT const> net);

public:

	/**
	 * @param net The network.
	 * @param seed The random number generator's seed.
	 */
	USBMPredictor(std::shared_ptr<NetworkT const> net, long int seed) :
			ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
					EdgesRandomOutputIteratorT>(net), seed(seed), rng(seed) {
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
	virtual double score(EdgeType const & e) {
		auto i = NetworkT::start(e);
		auto j = NetworkT::end(e);
		return scores[nodesMap[i]->num][nodesMap[j]->num];
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

}
/* namespace LinkPred */

#endif /* USBMPREDICTOR_HPP_ */
