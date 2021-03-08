/*
 * This file is part of LinkPred.
 *
 * LinkPred: A high performance library for link prediction in complex networks.
 * Copyright (C) 2021  by Said Kerrache.
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
 * @ingroup GraphAlg
 * @brief Contains the implementation of the optimization problem associated with matrix factorization.
 */

#ifndef MATFACTCG_HPP_
#define MATFACTCG_HPP_

#include "LinkPredConfig.hpp"
#include "linkpred/numerical/cg/cgdescent.hpp"
#include "linkpred/utils/randomgen.hpp"
#include <vector>

namespace LinkPred {

/**
 * @brief A simple structure to store matrix factoization problem data.
 */
struct MatFactPbData {
	long int i; /**< Start node. */
	long int j; /**< End node. */
	double target; /**< Target value (0/1 for unweighted graphs. */
};

/**
 * @brief Optimization problem associated with matrix factorization.
 */
class MatFactCG: public CG::CGDProblem {

	long int nbNodes; /**< The number of nodes. */
	std::vector<MatFactPbData> pbData; /**< The problem data. */
	int dim; /**< The dimension. */
	double lambda; /**< Regularization coeffcient. */
	double *coords; /**< The optimal coordinates found by solving MDS. */
	long int seed; /**< The seed. */
	RandomGen rng; /**< Random number generator. */
	int n; /**< Number of variables */
	int m; /**< Number of couples. */

protected:

#ifdef LINKPRED_WITH_OPENMP
	bool parallel = false; /**< To enable/disable parallelism. */
#endif

public:
	/**
	 * Constructor.
	 * @param nbNodes The number of nodes in the network.
	 * @param pbData The problem data.
	 * @param dim The problem dimensionality.
	 * @param lambda Regularization coeffcient.
	 * @param seed The random number generator's seed.
	 */
	MatFactCG(long int nbNodes, std::vector<MatFactPbData> const & pbData,
			int dim, double lambda, long int seed);

	/**
	 * Deleted copy constructor.
	 * @param that The object to copy.
	 */
	MatFactCG(const MatFactCG & that) = delete;

	/**
	 * Deleted copy assignment operator.
	 * @param that The object to copy.
	 */
	MatFactCG& operator=(const MatFactCG & that) = delete;

	/**
	 * Deleted move constructor.
	 * @param that The object to move.
	 */
	MatFactCG(MatFactCG && that) = delete;

	/**
	 * Deleted move assignment operator.
	 * @param that The object to move.
	 */
	MatFactCG& operator=(MatFactCG && that) = delete;

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
	virtual ~MatFactCG() = default;
};

} /* namespace LinkPred */

#endif /* MATFACTCG_HPP_ */
