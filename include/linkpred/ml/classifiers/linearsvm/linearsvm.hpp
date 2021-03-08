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
 * @brief Contains a wrapper of mlpack::smv::LinearSVM.
 */

#ifndef LINEARSVM_HPP_
#define LINEARSVM_HPP_

#include "LinkPredConfig.hpp"

#ifdef LINKPRED_WITH_MLPACK

#include "linkpred/ml/classifiers/classifier.hpp"
#include "linkpred/utils/randomgen.hpp"
#include <mlpack/core.hpp>
#include <mlpack/methods/linear_svm/linear_svm.hpp>

namespace LinkPred {


/**
 * @brief Linear SVM classifier.
 * @tparam InRndIt Input (features) iterator type.
 * @tparam OutRndIt Output (class) iterator type.
 * @tparam ScoreRndIt Classification scores iterator type.
 */
template<typename InRndIt = typename std::vector<Vec>::iterator,
		typename OutRndIt = typename std::vector<bool>::iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator> class LinearSVM: public mlpack::svm::LinearSVM<>,
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
	 * @param lambda Regularization coefficient.
	 * @param seed The random number generator's seed.
	 */
	LinearSVM() :
			Classifier<InRndIt, OutRndIt,
					ScoreRndIt>() {
		name = "LSVM";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	LinearSVM(LinearSVM const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	LinearSVM& operator =(LinearSVM const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	LinearSVM(LinearSVM &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	LinearSVM& operator =(LinearSVM &&that) = default;

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
	virtual ~LinearSVM() = default;

};


} /* namespace LinkPred */
#endif /* LINKPRED_WITH_MLPACK */
#endif /* LINEARSVM_HPP_ */
