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
 * @brief Contains the implementation of a link predictor that prestores the scores of edges.
 */

#ifndef DPSTPREDICTOR_HPP_
#define DPSTPREDICTOR_HPP_

#include <linkpred/predictors/directed/dlpredictor.hpp>
#include <memory>

namespace LinkPred {

/**
 * @brief A link predictor that prestores edge scores. This allows to seemingly integrate results from external link
 * prediction algorithms to LinkPred (for example, users may implement their own link prediction algorithm and then
 * use this link predictor to use compare their results to algorithms available in LinkPred).
 * @tparam Network The network type.
 * @tparam EdgeRndIt A random iterator type used to iterate on edges.
 * @tparam ScoreRndIt A random iterator type used to iterate on scores.
 */
template<typename Network = DNetwork<>,
		typename EdgeRndIt = typename std::vector<typename Network::Edge>::const_iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator,
		typename EdgeRndOutIt = typename std::vector<typename Network::Edge>::iterator> class DPSTPredictor: public DLPredictor<
		Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt> {

	using DLPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::net; /**< The network. */
#ifdef LINKPRED_WITH_OPENMP
	using DLPredictor<Network, EdgeRndIt, ScoreRndIt,
			EdgeRndOutIt>::parallel; /**< Whether the predictor runs in parallel. */
#endif
#ifdef LINKPRED_WITH_MPI
	using DLPredictor<Network, EdgeRndIt, ScoreRndIt,
			EdgeRndOutIt>::comm; /**< The MPI communicator. */
	using DLPredictor<Network, EdgeRndIt, ScoreRndIt,
	EdgeRndOutIt>::distributed; /**< Enable/disable distributed parallelism. */
#endif
	using DLPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::name; /**< The name of the predictor. */
	using NodeID = typename DLPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::NodeID; /**< The node IDs type. */
	using Edge = typename DLPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::Edge; /**< The edges type. */

protected:
	typename Network::template EdgeMapSP<double> edgeScores; /**< A map that stores edge scores. */

public:
	/**
	 * @param net The network.
	 */
	DPSTPredictor(std::shared_ptr<Network const> net) :
			DLPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>(net) {
		name = "PST";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	DPSTPredictor(DPSTPredictor const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	DPSTPredictor & operator =(DPSTPredictor const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	DPSTPredictor(DPSTPredictor && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	DPSTPredictor & operator =(DPSTPredictor && that) = default;

	/**
	 * Set the edge scores. This should contain the scores of all non-exisitng links of net.
	 * @param edgeScores A map contatining the scores of non-exisitng links of net.
	 */
	void setEdgeScores(
			typename Network::template EdgeMapSP<double> edgeScores) {
		this->edgeScores = edgeScores;
	}

	/**
	 * Load edge scores from file. This file should contain the scores of all non-exisitng links of net.
	 * @param fileName The name of the file containing edge scores. This should be a text file in which
	 * each line contains three columns separate by a space character (one or multiple spaces or tabs).
	 * The first two columns contain the the labels (not the internal IDs) of two nodes composing the edge.
	 * The third column contains the score. For example:
	 * A	B 	0.35
	 * or:
	 * 5	8	2.6
	 * All node labels that appear in this file must already exist in the network.
	 */
	void loadEdgeScores(std::string fileName);

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
	 * Compute the score of a single edge.
	 * @param e The edge.
	 * @return The score of e.
	 */
	virtual double score(Edge const & e);

	/**
	 * Destructor.
	 */
	virtual ~DPSTPredictor() = default;
};

} /* namespace LinkPred */

#endif /* DPSTPREDICTOR_HPP_ */
