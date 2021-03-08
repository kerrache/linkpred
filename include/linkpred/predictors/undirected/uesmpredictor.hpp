/*
 * This file is part of LinkPred.
 *
 * LinkPred: A high performance library for link prediction in complex networks.
 * Copyright (C) 2021  by Said Kerrache.
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
 * @brief Contains the implementation of an encoder-similarity measure link predictor.
 */

#ifndef UESMPREDICTOR_HPP_
#define UESMPREDICTOR_HPP_

#include <linkpred/predictors/undirected/ulpredictor.hpp>
#include <linkpred/graphalg/encoders/encoders.hpp>
#include <linkpred/ml/simmeasures/simmeasures.hpp>
#include <memory>

namespace LinkPred {

/**
 * @brief Encoder-Similarity measure link predictor.
 * @tparam Network The network type.
 * @tparam EdgeRndIt A random iterator type used to iterate on edges.
 * @tparam ScoreRndIt A random iterator type used to iterate on scores.
 */
template<typename Network = UNetwork<>,
		typename EdgeRndIt = typename std::vector<typename Network::Edge>::const_iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator,
		typename EdgeRndOutIt = typename std::vector<typename Network::Edge>::iterator> class UESMPredictor: public ULPredictor<
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

protected:
	std::shared_ptr<Encoder<Network>> encoder; /**< Network encoder. */
	std::shared_ptr<SimMeasure> simMeasure; /**< Similarity measure. */

public:
	/**
	 * @param net The network.
	 * @param encoder The encoder used to embed the network.
	 * @param simMeasure The similarity measure.
	 */
	UESMPredictor(std::shared_ptr<Network const> net,
			std::shared_ptr<Encoder<Network> > encoder,
			std::shared_ptr<SimMeasure> simMeasure) :
			ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>(net), encoder(
					encoder), simMeasure(simMeasure) {
		name = std::string("ESM") + "_" + encoder->getName() + "_"
				+ simMeasure->getName();
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	UESMPredictor(UESMPredictor const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	UESMPredictor& operator =(UESMPredictor const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	UESMPredictor(UESMPredictor &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	UESMPredictor& operator =(UESMPredictor &&that) = default;

	/**
	 * Initialize the predictor.
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
	virtual double score(Edge const &e);

	/**
	 * @return The encoder.
	 */
	const std::shared_ptr<Encoder<Network> >& getEncoder() const {
		return encoder;
	}

	/**
	 * Set the encoder.
	 * @param encoder The new encoder.
	 */
	void setEncoder(
			const std::shared_ptr<Encoder<Network> > &encoder) {
		this->encoder = encoder;
	}

	/**
	 * @return The similarity measure.
	 */
	const std::shared_ptr<SimMeasure>& getSimMeasure() const {
		return simMeasure;
	}

	/**
	 * Set the similarity measure.
	 * @param simMeasure The new similarity measure.
	 */
	void setSimMeasure(const std::shared_ptr<SimMeasure> &simMeasure) {
		this->simMeasure = simMeasure;
	}

	/**
	 * Destructor.
	 */
	virtual ~UESMPredictor() = default;

};
}
/* namespace LinkPred */

#endif /* UESMPREDICTOR_HPP_ */
