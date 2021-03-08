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
 * @brief Contains the interface of a classifier.
 */

#ifndef CLASSIFIER_HPP_
#define CLASSIFIER_HPP_

#include "LinkPredConfig.hpp"
#ifdef LINKPRED_WITH_OPENMP
#include <omp.h>
#endif
#ifdef LINKPRED_WITH_MPI
#include <mpi.h>
#endif
#include "linkpred/numerical/linear/vec.hpp"
#include <vector>

namespace LinkPred {

/**
 * @brief Interface of a binary classifier.
 * @tparam InRndIt Input (features) iterator type. Must be a random iterator to LinkPred::Vec.
 * @tparam OutRndIt Output (class) iterator type. Must be a random iterator to LinkPred::Vec.
 * @tparam ScoreRndIt Classification scores iterator type. Must be a random iterator to bool.
 */
template<typename InRndIt = typename std::vector<Vec>::iterator,
		typename OutRndIt = typename std::vector<bool>::iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator> class Classifier {

protected:
	std::string name; /**< The name of the classifier. */
#ifdef LINKPRED_WITH_OPENMP
	bool parallel = false; /**< Enable/disable shared-memory parallelism. */
#endif
#ifdef LINKPRED_WITH_MPI
	bool distributed = false; /**< Enable/disable distributed parallelism. */
	MPI_Comm comm = MPI_COMM_WORLD; /**< The MPI communicator. */
#endif

public:

	/**
	 * Default constructor.
	 */
	Classifier() = default;

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	Classifier(Classifier const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	Classifier& operator =(Classifier const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	Classifier(Classifier &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	Classifier& operator =(Classifier &&that) = default;

	/**
	 * Learn from data.
	 * @param trInBegin Iterator to the first example features (input).
	 * @param trInEnd Iterator to one-past-the-last example features (input).
	 * @param trOutBegin Iterator to the first example class (output).
	 * @param trOutEnd Iterator to one-past-the-last example class (output).
	 */
	virtual void learn(InRndIt trInBegin, InRndIt trInEnd, OutRndIt trOutBegin,
			OutRndIt trOutEnd) = 0;

	/**
	 * Predict.
	 * @param inBegin Iterator to the first instance features (input).
	 * @param inEnd Iterator to one-past-the-last instance features (input).
	 * @param scoresBegin Iterator to the first location where to store prediction scores. Memory must be pre-allocated.
	 */
	virtual void predict(InRndIt inBegin, InRndIt inEnd,
			ScoreRndIt scoresBegin) = 0;

	/**
	 * @return The name of the classifier.
	 */
	const std::string& getName() const {
		return name;
	}

	/**
	 * Set the name of the classifier.
	 * @param name The new name of the classifier.
	 */
	void setName(const std::string &name) {
		this->name = name;
	}

#ifdef LINKPRED_WITH_OPENMP
	/**
	 * @return Whether shared memory parallelism is enabled.
	 */
	bool isParallel() const {
		return parallel;
	}

	/**
	 * Enable/disable shared memory parallelism.
	 * @param parallel True to enable parallelism, false to disable it.
	 */
	void setParallel(bool parallel) {
		this->parallel = parallel;
	}
#endif

#ifdef LINKPRED_WITH_MPI
	/**
	 * @return Whether distributed memory parallelism is enabled.
	 */
	bool isDistributed() const {
		return distributed;
	}

	/**
	 * Enable/disable distributed memory parallelism.
	 * @param distributed True to enable distributed memory parallelism, false to disable it.
	 */
	void setDistributed(bool distributed) {
		this->distributed = distributed;
	}

	/**
	 * @return The MPI communicator.
	 */
	MPI_Comm getComm() const {
		return comm;
	}

	/**
	 * Set the MPI communicator.
	 * @param comm The new MPI communicator.
	 */
	void setComm(MPI_Comm const & comm) {
		this->comm = comm;
	}
#endif

	/**
	 * Destructor.
	 */
	virtual ~Classifier() = default;

};

} /* namespace LinkPred */

#endif /* CLASSIFIER_HPP_ */
