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
 * @brief Contains a wrapper of mlpack::maive_bayes::NaiveBayesClassifier.
 */

#ifndef NAIVEBAYES_HPP_
#define NAIVEBAYES_HPP_

#include "LinkPredConfig.hpp"

#ifdef LINKPRED_WITH_MLPACK

#include "linkpred/ml/classifiers/classifier.hpp"
#include "linkpred/utils/randomgen.hpp"
#include <mlpack/core.hpp>
#include <mlpack/methods/naive_bayes/naive_bayes_classifier.hpp>

namespace LinkPred {


/**
 * @brief Naive Bayes classifier.
 * @tparam InRndIt Input (features) iterator type.
 * @tparam OutRndIt Output (class) iterator type.
 * @tparam ScoreRndIt Classification scores iterator type.
 */
template<typename InRndIt = typename std::vector<Vec>::iterator,
		typename OutRndIt = typename std::vector<bool>::iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator> class NaiveBayes: public mlpack::naive_bayes::NaiveBayesClassifier<>,
		public Classifier<InRndIt, OutRndIt,
				ScoreRndIt> {

	using Classifier<InRndIt, OutRndIt,
			ScoreRndIt>::name; /**< Name of the classifier. */
#ifdef LINKPRED_WITH_OPENMP
	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::parallel; /**< Whether the predictor runs in parallel. */
#endif
#ifdef LINKPRED_WITH_MPI
	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::comm; /**< The MPI communicator. */
	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::distributed; /**< Enable/disable distributed parallelism. */
#endif

protected:
	std::size_t n = 0; /**< The number of observations (training examples).*/
	std::size_t m = 0; /**< The number of features. */

public:

	/**
	 * Constructor.
	 */
	NaiveBayes() :
			Classifier<InRndIt, OutRndIt,
					ScoreRndIt>() {
		name = "NVB";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	NaiveBayes(NaiveBayes const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	NaiveBayes& operator =(NaiveBayes const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	NaiveBayes(NaiveBayes &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	NaiveBayes& operator =(NaiveBayes &&that) = default;

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
	 * Destructor.
	 */
	virtual ~NaiveBayes() = default;

};


} /* namespace LinkPred */
#endif /* LINKPRED_WITH_MLPACK */
#endif /* NAIVEBAYES_HPP_ */
