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
 * @brief Contains the implementation of a solver of the multidimensional scaling problem
 * using the logarithmic (MULTISCALE) loss function.
 */

#ifndef MDS_HPP_
#define MDS_HPP_

#include <cstdlib>

namespace LinkPred {

/**
 * @brief MDS solution methods.
 */
enum MDSAlg {
	IpoptMDS, /**< Solution using Ipopt. */
	CGMDS /**< Solution using the CG Descent algorithm. */
};

/**
 * @brief Solve the MDS problem.
 */
class MDS {
protected:
	MDSAlg alg = CGMDS; /**< Algorithm for solving MDS */
	double tol = 1.0e-5; /**< Tolerance. */
	std::size_t nbRuns = 1; /**< Number of optimization runs. */

public:
	/**
	 * Constructor.
	 */
	MDS() = default;

	/**
	 * Constructor.
	 * @param alg The algorithm used to solve the MDS problem.
	 * @param tol The tolerance.
	 * @param nbRuns The number of times multidimensional scaling is run using different random initial positions.
	 */
	MDS(MDSAlg alg, double tol, std::size_t nbRuns) :
			alg(alg), tol(tol), nbRuns(nbRuns) {
	}

	/**
	 * @return The algorithm  used to solve the MDS problem.
	 */
	MDSAlg getAlg() {
		return alg;
	}

	/**
	 * Set the MDS algorithm.
	 * @param alg The algorithm  used to solve the MDS problem.
	 */
	void setAlg(MDSAlg alg) {
		MDS::alg = alg;
	}

	/**
	 * @return The tolerance.
	 */
	double getTol() {
		return tol;
	}

	/**
	 * Set the tolerance.
	 * @param tol The new value of the tolerance.
	 */
	void setTol(double tol) {
		MDS::tol = tol;
	}

	/**
	 * @return The number of times multidimensional scaling is run using different random initial positions.
	 */
	std::size_t getNbRuns() {
		return nbRuns;
	}

	/**
	 * Set the number of runs.
	 * @param nbRuns The new value of the number of runs.
	 */
	void setNbRuns(std::size_t nbRuns) {
		MDS::nbRuns = nbRuns;
	}

	/**
	 * Solve the MDS problem.
	 * @param nbNodes Number of nodes (or more generally, the number of data points).
	 * @param nbKnownCouples Number of couples for which the target distance is specified.
	 * @param sqDist Upper triangular matrix of squared distances without the diagonal stored sequentially as one vector.
	 * @param weight The weight associated with every distance (can be zero, if the distance is not relevant or unknown).
	 * @param dim The problem dimensionality.
	 * @param coord The coordinates (this is the solution to the problem).
	 * @param seed The random number generator's seed.
	 * @param init Whether to initialize points randomly or use the current values of the coordinates as the initial positions of the data points.
	 * @return The objective function value.
	 */
	double solve(std::size_t nbNodes, std::size_t nbKnownCouples,
			double *sqDist, double *weight, std::size_t dim, double *coord,
			long int seed, bool init = true);

#ifdef WITH_IPOPT
	/**
	 * MDS using Ipopt.
	 * @param nbNodes Number of nodes (or more generally, the number of data points).
	 * @param nbKnownCouples Number of couples for which the target distance is specified.
	 * @param sqDist Upper triangular matrix of squared distances without the diagonal stored sequentially as one vector.
	 * @param weight The weight associated with every distance (can be zero, if the distance is not relevant or unknown).
	 * @param dim The problem dimensionality.
	 * @param coords The coordinates (this is the solution to the problem).
	 * @param seed The random number generator's seed.
	 * @param init Whether to initialize points randomly or use the current values of the coordinates as the initial positions of the data points.
	 * @return The objective function value.
	 */
	double ipoptMDS(std::size_t nbNodes, std::size_t nbKnownCouples,
			std::size_t dim, double *sqDist, double *weight, double *coord,
			bool init);
#endif
	/**
	 * MDS using CG.
	 * @param nbNodes Number of nodes (or more generally, the number of data points).
	 * @param nbKnownCouples Number of couples for which the target distance is specified.
	 * @param sqDist Upper triangular matrix of squared distances without the diagonal stored sequentially as one vector.
	 * @param weight The weight associated with every distance (can be zero, if the distance is not relevant or unknown).
	 * @param dim The problem dimensionality.
	 * @param coord The coordinates (this is the solution to the problem).
	 * @param seed The random number generator's seed.
	 * @param init Whether to initialize points randomly or use the current values of the coordinates as the initial positions of the data points.
	 * @return The objective function value.
	 */
	double cgMDS(std::size_t nbNodes, std::size_t nbKnownCouples,
			std::size_t dim, double *sqDist, double *weight, double *coord,
			long int seed, bool init);

	/**
	 * Destructor.
	 */
	virtual ~MDS() = default;
};

} /* namespace LinkPred */

#endif /* INCLUDE_MDS_HPP_ */
