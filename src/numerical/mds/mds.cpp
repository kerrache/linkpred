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

#include "linkpred/numerical/mds/mds.hpp"
#ifdef WITH_IPOPT
#include "linkpred/IpIpoptApplication.hpp"
#include "linkpred/IpSolveStatistics.hpp"
#include "linkpred/logmdsnlp.hpp"
#endif
#include "linkpred/numerical/mds/logmdscg.hpp"
#include "linkpred/numerical/cg/cgdescent.hpp"
#include "linkpred/utils/log.hpp"
#include <linkpred/utils/miscutils.hpp>
#include <limits>

namespace LinkPred {

double MDS::solve(std::size_t nbNodes, std::size_t nbKnownCouples,
		double *sqDist, double *weight, std::size_t dim, double *coords,
		long int seed, bool init) {
	switch (alg) {
#ifdef WITH_IPOPT
	case IpoptMDS:
	return ipoptMDS(nbNodes, nbKnownCouples, dim, sqDist, weight, coords,
			init);
	break;
#endif
	case CGMDS:
		return cgMDS(nbNodes, nbKnownCouples, dim, sqDist, weight, coords, seed,
				init);
		break;

	default:
		logger(logError, "Unknown MDS solution method. Abort.")
		exit(1);
	}
	return std::numeric_limits<double>::max();
}

#ifdef WITH_IPOPT
double MDS::ipoptMDS(std::size_t nbNodes, std::size_t nbKnownCouples, std::size_t dim, double *sqDist,
		double *weight, double *coords, bool init) {
	logger(logError, "Performing MDS using Ipopt...")
	double bestObj = std::numeric_limits<double>::max();
	std::size_t n = dim * nbNodes;
	double* bestCoords = new double[n];
	bool copyBack = false;

	for (std::size_t i = 0; i < nbRuns; i++) {
		auto mdsPb = new LogMDSNlp(nbNodes, nbKnownCouples, sqDist, weight, dim,
				coords, init);
		Ipopt::SmartPtr<TNLP> mynlp = mdsPb;

		Ipopt::SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
		app->Options()->SetNumericValue("tol", tol);
		app->Options()->SetStringValue("linear_solver", "MA57");
		app->Options()->SetStringValue("hessian_approximation",
				"limited-memory");

		// Initialize the IpoptApplication and process the options
		Ipopt::ApplicationReturnStatus status;
		status = app->Initialize();
		if (status != Solve_Succeeded) {
			logger(logError, "Error during Ipopt initialization!")
			return std::numeric_limits<double>::max();
		}

		// Solve
		status = app->OptimizeTNLP(mynlp);

		if (status != Solve_Succeeded) {
			logger(logInfo, "Ipopt cannot solve the problem!")
		}

		double obj = mdsPb->getObj();
		std::cout << "MDS Obj: " << obj << std::endl;
		if (obj < bestObj) {
			bestObj = obj;
			for (std::size_t i = 0; i < n; i++) {
				bestCoords[i] = coords[i];
			}
			copyBack = false;
		} else {
			copyBack = true;
		}
	}

	if (copyBack) {
		for (std::size_t i = 0; i < n; i++) {
			coords[i] = bestCoords[i];
		}
	}

	delete[] bestCoords;
	logger(logDebug, "Done")
	return bestObj;
}
#endif

double MDS::cgMDS(std::size_t nbNodes, std::size_t nbKnownCouples,
		std::size_t dim, double *sqDist, double *weight, double *coords,
		long int seed, bool init) {
	logger(logDebug, "Performing MDS using CG...")
	double bestObj = std::numeric_limits<double>::max();
	std::size_t n = dim * nbNodes;
	double* bestCoords = new double[n];
	bool copyBack = false;
	CG::cg_stats stats;

	for (std::size_t i = 0; i < nbRuns; i++) {
		CG::CGDProblem* pb = new LogMDSCG(nbNodes, nbKnownCouples, sqDist,
				weight, dim, coords, seed);

		CG::CGDescent solver(pb);
		solver.cg_descent(coords, n, &stats, nullptr, tol, nullptr, init);

		double obj = stats.f;
		if (obj < bestObj) {
			bestObj = obj;
			for (std::size_t i = 0; i < n; i++) {
				bestCoords[i] = coords[i];
			}
			copyBack = false;
		} else {
			copyBack = true;
		}
		delete pb;
	}

	if (copyBack) {
		for (std::size_t i = 0; i < n; i++) {
			coords[i] = bestCoords[i];
		}
	}

	delete[] bestCoords;
	logger(logDebug, "Done")
	return bestObj;
}

} /* namespace LinkPred */
