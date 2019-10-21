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
 * @brief Contains the implementation of a logistic regression algorithm.
 */

#ifndef LOGISTICREGRESSER_HPP_
#define LOGISTICREGRESSER_HPP_

#include "linkpred/utils/randomgen.hpp"
#include "linkpred/numerical/linear/vec.hpp"
#include <vector>

namespace LinkPred {

/**
 * @brief Logistic regression algorithm.
 * @tparam InputRandomIterator Input (features) iterator type.
 * @tparam OutputRandomIterator Output (class) iterator type.
 * @tparam ScoreRandomIterator Classification scores iterator type.
 */
template<typename InputRandomIterator = typename std::vector<Vec>::iterator,
		typename OutputRandomIterator = typename std::vector<int>::iterator,
		typename ScoreRandomIterator = typename std::vector<double>::iterator> class LogisticRegresser {

	double lambda; /**< Regularization coefficient. */
	long int seed; /**< Seed of the random number generator. */
	RandomGen rng; /**< A random number generator. */
	std::size_t n = 0; /**< The number of observations (training examples).*/
	std::size_t m = 0; /**< The number of features. */
	double tol = 1.0e-6; /**< Tolerance used to stop the optimization procedure. */
	Vec theta; /**< The model parameters. */
	std::vector<double> maxF; /**< Stores the maximum value of every feature (used to scale features). */
	std::vector<double> minF; /**< Stores the minimum value of every feature (used to scale features). */
	bool parallel = false; /**< Enable/disable parallelism. */

public:

	/**
	 * Constructor.
	 * @param lambda Regularization coefficient.
	 * @param seed The random number generator's seed.
	 */
	LogisticRegresser(double lambda, long int seed) :
			lambda(lambda), seed(seed), rng(seed) {
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	LogisticRegresser(LogisticRegresser const & that);

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	LogisticRegresser & operator =(LogisticRegresser const & that);

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	LogisticRegresser(LogisticRegresser && that);

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	LogisticRegresser & operator =(LogisticRegresser && that);

	/**
	 * Learn from data.
	 * @param trainingInputsBegin Iterator to the first example features (input).
	 * @param trainingInputsEnd Iterator to one-past-the-last example features (input).
	 * @param trainingOutputsBegin Iterator to the first example class (output).
	 * @param trainingOutputsEnd Iterator to one-past-the-last example class (output).
	 */
	void learn(InputRandomIterator trainingInputsBegin,
			InputRandomIterator trainingInputsEnd,
			OutputRandomIterator trainingOutputsBegin,
			OutputRandomIterator trainingOutputsEnd);

	/**
	 * Predict.
	 * @param testInputsBegin Iterator to the first instance features (input).
	 * @param testInputsEnd Iterator to one-past-the-last instance features (input).
	 * @param scoresBegin Iterator to the first location where to store prediction scores.
	 */
	void predict(InputRandomIterator testInputsBegin,
			InputRandomIterator testInputsEnd,
			ScoreRandomIterator scoresBegin) const;

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
	 * @return The model paramneters.
	 */
	const Vec& getTheta() const {
		return theta;
	}

	/**
	 * @return Whether parallelism is enabled.
	 */
	bool isParallel() const {
		return parallel;
	}

	/**
	 * Enable/disable parallelism.
	 * @param parallel True to enable parallelism, false to disable it.
	 */
	void setParallel(bool parallel) {
		this->parallel = parallel;
	}

	/**
	 * Destructor.
	 */
	virtual ~LogisticRegresser() = default;

};

} /* namespace LinkPred */

#endif /* LOGISTICREGRESSER_HPP_ */
