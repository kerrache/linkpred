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
 * @brief Contains the implementation of the optimization problem of computing w.
 */

#ifndef INCLUDE_KABLAMBDACG_HPP_
#define INCLUDE_KABLAMBDACG_HPP_

#include "LinkPredConfig.hpp"
#include "linkpred/numerical/cg/cgdescent.hpp"
#ifdef WITH_OPENMP
#include <omp.h>
#endif
#include <vector>
#include <stdexcept>

namespace LinkPred {

/**
 * @brief Finds best w using conjugate gradient.
 */
class KABLambdaCG: public CG::CGDProblem {
protected:
#ifdef WITH_OPENMP
	bool parallel = false; /**< Run in parallel ? */
#endif
	double sol = 0.0; /**< The value of w computed. */
	std::vector<double> posSUM; /**< Positive scores. */
	std::vector<double> negSUM; /**< Negative scores. */
	std::vector<double> posNSC; /**< Positive scores. */
	std::vector<double> negNSC; /**< Negative scores. */
	std::size_t nbPos; /**< Number of positive links. */
	std::size_t nbNeg; /**< Number of negative links. */
	std::size_t nb; /**< Total: nbPos + nbNeg. */
	double w; /**< Weight. */

public:

	/**
	 * Default constructor.
	 */
	KABLambdaCG();

	/**
	 * Constructor.
	 */
	KABLambdaCG(std::vector<double> const & posSUM,
			std::vector<double> const & negSUM,
			std::vector<double> const & posNSC,
			std::vector<double> const & negNSC, double w) :
			posSUM(posSUM), negSUM(negSUM), posNSC(posNSC), negNSC(negNSC), w(w) {
		nbPos = posSUM.size();
		if (nbPos != posNSC.size()) {
			throw std::invalid_argument("Vectors must have the same size");
		}
		nbNeg = negSUM.size();
		if (nbNeg != negNSC.size()) {
			throw std::invalid_argument("Vectors must have the same size");
		}
		nb = nbPos + nbNeg;
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	KABLambdaCG(KABLambdaCG const & that) = delete;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	KABLambdaCG & operator =(KABLambdaCG const & that) = delete;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	KABLambdaCG(KABLambdaCG && that) = delete;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	KABLambdaCG & operator =(KABLambdaCG && that) = delete;

#ifdef WITH_OPENMP
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
#endif

	/**
	 * Initializes x.
	 * @param x The variables.
	 * @param n The size.
	 */
	virtual void init(double * x, CG::INT n);

	/**
	 * Objective function.
	 * @param x The variables.
	 * @param n The size.
	 * @return The value of the objective.
	 */
	virtual double f(double * x, CG::INT n);

	/**
	 * Gradient.
	 * @param grad_f The gradient (output).
	 * @param x The variables.
	 * @param n The size.
	 */
	virtual void grad(double * grad_f, double * x, CG::INT n);

	/**
	 * Compute objective and gradient at the same time. This is a default
	 * implementation that should be overloaded if needed.
	 * @param grad_f The gradient (output).
	 * @param x The variables.
	 * @param n The size.
	 * @return The value of the objective.
	 */
	virtual double fgrad(double * grad_f, double * x, CG::INT n);

	/**
	 * Finalize the solution.
	 * @param x The variables.
	 * @param n The size.
	 */
	virtual void finalize(double * x, CG::INT n);

	/**
	 * @return The solution.
	 */
	double getSol() const {
		return sol;
	}

	/**
	 * Destructor.
	 */
	virtual ~KABLambdaCG() = default;
};

} /* namespace LinkPred */

#endif /* INCLUDE_KABLAMBDACG_HPP_ */
