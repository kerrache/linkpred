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
 * @brief Contains the implementation of a logistic regression algorithm.
 */

#ifndef LOGISTICREGRESSER_HPP_
#define LOGISTICREGRESSER_HPP_

#include "linkpred/ml/classifiers/classifier.hpp"
#include "linkpred/utils/randomgen.hpp"

namespace LinkPred {

/**
 * @brief Logistic regression algorithm.
 * @tparam InRndIt Input (features) iterator type.
 * @tparam OutRndIt Output (class) iterator type.
 * @tparam ScoreRndIt Classification scores iterator type.
 */
template<typename InRndIt = typename std::vector<Vec>::iterator,
		typename OutRndIt = typename std::vector<bool>::iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator> class LogisticRegresser: public Classifier<
		InRndIt, OutRndIt, ScoreRndIt> {

	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::name; /**< Name of the classifier. */
#ifdef LINKPRED_WITH_OPENMP
	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::parallel; /**< Whether the predictor runs in parallel. */
#endif
#ifdef LINKPRED_WITH_MPI
	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::comm; /**< The MPI communicator. */
	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::distributed; /**< Enable/disable distributed parallelism. */
#endif

	double lambda; /**< Regularization coefficient. */
	long int seed; /**< Seed of the random number generator. */
	RandomGen rng; /**< A random number generator. */
	std::size_t n = 0; /**< The number of observations (training examples).*/
	std::size_t m = 0; /**< The number of features. */
	double tol = 1.0e-6; /**< Tolerance used to stop the optimization procedure. */
	Vec theta; /**< The model parameters. */
	std::vector<double> maxF; /**< Stores the maximum value of every feature (used to scale features). */
	std::vector<double> minF; /**< Stores the minimum value of every feature (used to scale features). */

	/**
	 * Stores range in a vector and adds the constant feature 1.
	 */
	std::vector<Vec> toOneVec(InRndIt begin, InRndIt end);

public:

	/**
	 * Constructor.
	 * @param lambda Regularization coefficient.
	 * @param seed The random number generator's seed.
	 */
	LogisticRegresser(double lambda, long int seed) :
			Classifier<InRndIt, OutRndIt, ScoreRndIt>(), lambda(lambda), seed(
					seed), rng(seed) {
		name = "LGR";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	LogisticRegresser(LogisticRegresser const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	LogisticRegresser& operator =(LogisticRegresser const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	LogisticRegresser(LogisticRegresser &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	LogisticRegresser& operator =(LogisticRegresser &&that) = default;

	/**
	 * Learn from data.
	 * @param trInBegin Iterator to the first example features (input).
	 * @param trInEnd Iterator to one-past-the-last example features (input).
	 * @param trOutBegin Iterator to the first example class (output).
	 * @param trOutEnd Iterator to one-past-the-last example class (output).
	 */
	virtual void learn(InRndIt trInBegin, InRndIt trInEnd, OutRndIt trOutBegin,
			OutRndIt trOutEnd);

	/**
	 * Predict.
	 * @param inBegin Iterator to the first instance features (input).
	 * @param inEnd Iterator to one-past-the-last instance features (input).
	 * @param scoresBegin Iterator to the first location where to store prediction scores.
	 */
	virtual void predict(InRndIt inBegin, InRndIt inEnd,
			ScoreRndIt scoresBegin);

	/**
	 * @return The regularization coefficient.
	 */
	double getLambda() const {
		return lambda;
	}

	/**
	 * Set the regularization coefficient.
	 * @param lambda The new regularization coefficient.
	 */
	void setLambda(double lambda) {
		this->lambda = lambda;
	}

	/**
	 * @return The tolerance.
	 */
	double getTol() const {
		return tol;
	}

	/**
	 * Set the tolerance.
	 * @param tol The new value of the tolerance.
	 */
	void setTol(double tol) {
		this->tol = tol;
	}

	/**
	 * @return The model parameters.
	 */
	const Vec& getTheta() const {
		return theta;
	}

	/**
	 * Destructor.
	 */
	virtual ~LogisticRegresser() = default;

};

} /* namespace LinkPred */

#endif /* LOGISTICREGRESSER_HPP_ */
