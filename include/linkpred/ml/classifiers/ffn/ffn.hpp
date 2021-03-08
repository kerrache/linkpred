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
 * @ingroup ML
 * @brief Contains a wrapper of mlpack::ann::FFN (feed-forward neural network).
 */

#ifndef FFN_HPP_
#define FFN_HPP_

#include "LinkPredConfig.hpp"

#ifdef LINKPRED_WITH_MLPACK

#include "linkpred/ml/classifiers/classifier.hpp"
#include "linkpred/utils/randomgen.hpp"
#include <mlpack/core.hpp>
#include <mlpack/methods/ann/layer/layer.hpp>
#include <mlpack/methods/ann/ffn.hpp>
#include <cmath>

namespace LinkPred {


/**
 * @brief Feed-forward neural network.
 * @tparam InRndIt Input (features) iterator type.
 * @tparam OutRndIt Output (class) iterator type.
 * @tparam ScoreRndIt Classification scores iterator type.
 */
template<typename InRndIt = typename std::vector<Vec>::iterator,
		typename OutRndIt = typename std::vector<bool>::iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator> class FFN: public mlpack::ann::FFN<>,
		public Classifier<InRndIt, OutRndIt,
				ScoreRndIt> {

	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::name; /**< Name of the classifier. */
#ifdef LINKPRED_WITH_OPENMP
	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::parallel; /**< Whether the predictor runs in parallel. */
#endif
#ifdef LINKPRED_WITH_MPI
	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::comm; /**< The MPI communicator. */
	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::distributed; /**< Enable/disable distributed parallelism. */
#endif

private:
	std::size_t n = 0; /**< The number of observations (training examples).*/
	std::size_t m = 0; /**< The number of features. */

	/**
	 * Sigmoid
	 * @param x input.
	 * @return sigmoid of x: 1/(1+exp(-x)).
	 */
	inline double sigmoid(double x) {
		return 1.0 / (1.0 + std::exp(-x));
	}

public:

	/**
	 * Constructor.
	 * @param lambda Regularization coefficient.
	 * @param seed The random number generator's seed.
	 */
	FFN() :	Classifier<InRndIt, OutRndIt,
					ScoreRndIt>() {
		name = "FFN";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	FFN(FFN const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	FFN& operator =(FFN const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	FFN(FFN &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	FFN& operator =(FFN &&that) = default;

	/**
	 * Learn from data.
	 * @param trInBegin Iterator to the first example features (input).
	 * @param trInEnd Iterator to one-past-the-last example features (input).
	 * @param trOutBegin Iterator to the first example class (output).
	 * @param trOutEnd Iterator to one-past-the-last example class (output).
	 */
	virtual void learn(InRndIt trInBegin,
			InRndIt trInEnd,
			OutRndIt trOutBegin,
			OutRndIt trOutEnd);

	/**
	 * Predict.
	 * @param inBegin Iterator to the first instance features (input).
	 * @param inEnd Iterator to one-past-the-last instance features (input).
	 * @param scoresBegin Iterator to the first location where to store prediction scores.
	 */
	virtual void predict(InRndIt inBegin,
			InRndIt inEnd, ScoreRndIt scoresBegin);

	/**
	 * Set automatically the architecture.
	 * @param dim Input dimension (number of features).
	 */
	void setAutoArch(int dim);

	/**
	 * Destructor.
	 */
	virtual ~FFN() = default;

};


} /* namespace LinkPred */
#endif /* LINKPRED_WITH_MLPACK */
#endif /* FFN_HPP_ */
