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
 * @brief Contains the implementation of an encoder-classifier link predictor.
 */

#ifndef UECLPREDICTOR_HPP_
#define UECLPREDICTOR_HPP_

#include <linkpred/predictors/undirected/ulpredictor.hpp>
#include <linkpred/graphalg/encoders/encoders.hpp>
#include <linkpred/ml/classifiers/classifiers.hpp>
#include <memory>

namespace LinkPred {

/**
 * @brief Encoder-classifier link predictor.
 * @tparam Network The network type.
 * @tparam EdgeRndIt A random iterator type used to iterate on edges.
 * @tparam ScoreRndIt A random iterator type used to iterate on scores.
 */
template<typename Network = UNetwork<>,
		typename EdgeRndIt = typename std::vector<typename Network::Edge>::const_iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator,
		typename EdgeRndOutIt = typename std::vector<typename Network::Edge>::iterator> class UECLPredictor: public ULPredictor<
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
	std::shared_ptr<Classifier<>> classifier; /**< Classifier. */
	double posRatio = 1.0; /**< Ratio of positive edges used in the training of the classifier. */
	double negRatio = 1.0; /**< Ratio of negative edges used in the training of the classifier. */

private:
	RandomGen rng; /**< Random number generator. */

public:
	/**
	 * @param net The network.
	 * @param encoder The encoder used to embed the network.
	 * @param classifier The classifier used to discrminate between positive and negative links.
	 * @param seed Seed for the random number generator.
	 */
	UECLPredictor(std::shared_ptr<Network const> net,
			std::shared_ptr<Encoder<Network> > encoder,
			std::shared_ptr<Classifier<> > classifier, long int seed) :
			ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>(net), encoder(
					encoder), classifier(classifier), rng(seed) {
		name = std::string("ENC") + "_" + encoder->getName() + "_"
				+ classifier->getName();
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	UECLPredictor(UECLPredictor const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	UECLPredictor& operator =(UECLPredictor const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	UECLPredictor(UECLPredictor &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	UECLPredictor& operator =(UECLPredictor &&that) = default;

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
	 * @return The classifier.
	 */
	const std::shared_ptr<Classifier<> >& getClassifier() const {
		return classifier;
	}

	/**
	 * @return The encoder.
	 */
	const std::shared_ptr<Encoder<Network> >& getEncoder() const {
		return encoder;
	}

	/**
	 * Set the classifier.
	 * @param classifier The new classifier.
	 */
	void setClassifier(const std::shared_ptr<Classifier<> > &classifier) {
		this->classifier = classifier;
	}

	/**
	 * Set the encoder.
	 * @param encoder The new encoder.
	 */
	void setEncoder(const std::shared_ptr<Encoder<Network> > &encoder) {
		this->encoder = encoder;
	}

	/**
	 * @return Ratio of negative edges used in the training of the classifier.
	 */
	double getNegRatio() const {
		return negRatio;
	}

	/**
	 * Set the ratio of negative edges used in the training of the classifier.
	 * @param negRatio Ratio of negative edges used in the training of the classifier.
	 */
	void setNegRatio(double negRatio) {
		this->negRatio = negRatio;
	}

	/**
	 * @return Ratio of positive edges used in the training of the classifier.
	 */
	double getPosRatio() const {
		return posRatio;
	}

	/**
	 * Set the ratio of positive edges used in the training of the classifier.
	 * @param posRatio Ratio of positive edges used in the training of the classifier.
	 */
	void setPosRatio(double posRatio) {
		this->posRatio = posRatio;
	}

	/**
	 * Destructor.
	 */
	virtual ~UECLPredictor() = default;

};
}
/* namespace LinkPred */

#endif /* UECLPREDICTOR_HPP_ */
