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
 * @brief Contains a random classiffier (for debugging purposes).
 */

#ifndef RNDCLASSIFIER_HPP_
#define RNDCLASSIFIER_HPP_

#include "LinkPredConfig.hpp"

#include "linkpred/ml/classifiers/classifier.hpp"
#include "linkpred/utils/randomgen.hpp"

namespace LinkPred {

/**
 * @brief Random classifier.
 * @tparam InRndIt Input (features) iterator type.
 * @tparam OutRndIt Output (class) iterator type.
 * @tparam ScoreRndIt Classification scores iterator type.
 */
template<typename InRndIt = typename std::vector<Vec>::iterator,
		typename OutRndIt = typename std::vector<bool>::iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator> class RndClassifier: public Classifier<
		InRndIt, OutRndIt, ScoreRndIt> {

	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::name; /**< Name of the classifier. */
#ifdef LINKPRED_WITH_OPENMP
	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::parallel; /**< Whether the predictor runs in parallel. */
#endif
#ifdef LINKPRED_WITH_MPI
	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::comm; /**< The MPI communicator. */
	using Classifier<InRndIt, OutRndIt, ScoreRndIt>::distributed; /**< Enable/disable distributed parallelism. */
#endif

protected:
	RandomGen rng; /**< Random number generator.*/

public:

	/**
	 * Constructor.
	 * @param seed The random number generator's seed.
	 */
	RndClassifier(long int seed) :
			Classifier<InRndIt, OutRndIt, ScoreRndIt>(), rng(seed) {
		name = "RND";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	RndClassifier(RndClassifier const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	RndClassifier& operator =(RndClassifier const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	RndClassifier(RndClassifier &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	RndClassifier& operator =(RndClassifier &&that) = default;

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
	 * Destructor.
	 */
	virtual ~RndClassifier() = default;

};

} /* namespace LinkPred */
#endif /* RNDCLASSIFIER_HPP_ */
