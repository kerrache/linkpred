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
 * @ingroup Numerical
 * @brief Contains the implementation of the optimization problem associated with multidimensional scaling
 * using the logarithmic (MULTISCALE) loss function.
 */

#ifndef LOGMDSCG_HPP_
#define LOGMDSCG_HPP_

#include "linkpred/numerical/cg/cgdescent.hpp"
#include "linkpred/utils/randomgen.hpp"

namespace LinkPred {

/**
 * @brief Optimization problem associated with multidimensional scaling
 * using the logarithmic (MULTISCALE) loss function
 */
class LogMDSCG: public CG::CGDProblem {

	std::size_t nbNodes; /**< Number of nodes. */
	std::size_t nbKnownCouples; /**< Number of known couples. */
	double *sqDist; /**< The squared distances. */
	double *weight; /**< The weights. */
	std::size_t dim; /**< The dimension. */
	double *coords; /**< The optimal coordinates found by solving MDS. */
	std::size_t n; /**< Number of variables */
	double scl = 1; /**< Scaling factor. */
	long int seed; /**< The seed. */
	RandomGen rng; /**< Random number generator. */

	/**
	 * Scale the problem.
	 */
	void scale();

protected:
	bool parallel = false; /**< To enable/disable parallelism. */

public:
	/**
	 * Constructor.
	 * @param nbNodes Number of nodes (or more generally, the number of data points).
	 * @param nbKnownCouples Number of couples for which the target distance is specified.
	 * @param sqDist Upper triangular matrix of squared distances without the diagonal stored sequentially as one vector.
	 * @param weight The weight associated with every distance (can be zero, if the distance is not relevant or unknown).
	 * @param dim The problem dimensionality.
	 * @param coords The coordinates (this is the solution to the problem).
	 * @param seed The random number generator's seed.
	 */
	LogMDSCG(std::size_t nbNodes, std::size_t nbKnownCouples, double * sqDist,
			double * weight, std::size_t dim, double * coords, long int seed);

	/**
	 * Deleted copy constructor.
	 * @param that The object to copy.
	 */
	LogMDSCG(const LogMDSCG & that) = delete;

	/**
	 * Deleted copy assignment operator.
	 * @param that The object to copy.
	 */
	LogMDSCG& operator=(const LogMDSCG & that) = delete;

	/**
	 * Deleted move constructor.
	 * @param that The object to move.
	 */
	LogMDSCG(LogMDSCG && that) = delete;

	/**
	 * Deleted move assignment operator.
	 * @param that The object to move.
	 */
	LogMDSCG& operator=(LogMDSCG && that) = delete;

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
	 * Default destructor.
	 */
	virtual ~LogMDSCG() = default;
};

} /* namespace LinkPred */

#endif /* LOGMDSCG_HPP_ */
