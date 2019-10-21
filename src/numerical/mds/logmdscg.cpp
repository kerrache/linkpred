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

#include "linkpred/numerical/mds/logmdscg.hpp"
#include "linkpred/utils/randomgen.hpp"
#include <iostream>
#include <omp.h>
#include "linkpred/utils/utilities.hpp"

// TODO flatten for loops and add omp parallelization
namespace LinkPred {

LogMDSCG::LogMDSCG(std::size_t nbNodes, std::size_t nbKnownCouples,
		double * sqDist, double * weight, std::size_t dim, double * coords,
		long int seed) :
		seed(seed), rng(seed) {

	this->nbNodes = nbNodes;
	this->nbKnownCouples = nbKnownCouples;
	this->sqDist = sqDist;
	this->weight = weight;
	this->dim = dim;
	this->coords = coords;
	n = dim * nbNodes; // Number of variables
	scale();
}

void LogMDSCG::scale() {

	double meanD = 0;
	double sumW = 0;
	//The following for loop is not parallelized, because it is too sensitive. Error will propagate.
////#pragma omp parallel for reduction(+:meanD), reduction(+:sumW)
	for (std::size_t i = 0; i < nbNodes * (nbNodes - 1) / 2; i++) {
		if (weight[i] > 0) {
			meanD += weight[i] * std::log(sqDist[i]);
			sumW += weight[i];
		}
	}

	scl = std::exp(meanD / sumW);
//#pragma omp parallel for
	for (std::size_t i = 0; i < nbNodes * (nbNodes - 1) / 2; i++) {
		sqDist[i] /= scl;
	}

	double sqScl = std::sqrt(scl);
//#pragma omp parallel for
	for (std::size_t i = 0; i < nbNodes * dim; i++) {
		coords[i] /= sqScl;
	}

}

void LogMDSCG::init(double * x, CG::INT n) {
	double maxD = 0;
//#pragma omp parallel for reduction(max:maxD)
	for (std::size_t i = 0; i < nbNodes * (nbNodes - 1) / 2; i++) {
		if (sqDist[i] > maxD) {
			maxD = sqDist[i];
		}
	}
	maxD = std::sqrt(maxD);

//////#pragma omp parallel for // TODO Check if safe
	RandomGen rng2;
	double lambda = 1.0e-10;
	for (long int i = 0; i < n; i++) {
		x[i] = rng.getDouble(0, maxD) + lambda * rng2.getDouble(0, 1);
	}
}

double LogMDSCG::f(double * x, CG::INT n) {

	std::size_t l = 0;
	double d = 0.0;
	double sum = 0.0;
	for (std::size_t i = 0; i < nbNodes; i++) {
		for (std::size_t j = i + 1; j < nbNodes; j++, l++) {
			double w = weight[l];
			if (w > 0) { //i.e. known couple
				double ksi = sqDist[l]; // Observed distance
				d = 0.0;
				for (std::size_t k = 0; k < dim; k++) { // Compute the current distance between the two nodes
					double v = x[i * dim + k] - x[j * dim + k];
					d += v * v;
				}
				double v = std::log(d / ksi);
				sum += w * v * v;
			}
		}
	}

	return (sum / nbKnownCouples);
}

void LogMDSCG::grad(double * grad_f, double * x, CG::INT n) {
// Initialize to 0
//#pragma omp parallel for
	for (long int i = 0; i < n; i++) {
		grad_f[i] = 0.0;
	}

// Go over all couples
	std::size_t l = 0;
	for (std::size_t i = 0; i < nbNodes; i++) {
		for (std::size_t j = i + 1; j < nbNodes; j++, l++) {
			double ksi = sqDist[l]; // Observed distance
			double w = weight[l];
//			std::cout << w << "\t";
			if (w > 0) { //i.e. known couple
				double d = 0.0;
				for (std::size_t k = 0; k < dim; k++) { // Compute the current distance between the two nodes
					std::size_t indexI = i * dim + k;
					std::size_t indexJ = j * dim + k;
					double v = x[indexI] - x[indexJ];
					d += v * v;
				}
				double val = std::log(d / ksi);
				for (std::size_t k = 0; k < dim; k++) {
					std::size_t indexI = i * dim + k;
					std::size_t indexJ = j * dim + k;

					double v = w * (2.0 * val) * (x[indexI] - x[indexJ])
							/ (d * nbKnownCouples);
					grad_f[indexI] += v;
					grad_f[indexJ] -= v;
				}
			}
		}
	}
}

double LogMDSCG::fgrad(double * grad_f, double * x, CG::INT n) {
// Initialize to 0
//#pragma omp parallel for
	for (long int i = 0; i < n; i++) {
		grad_f[i] = 0.0;
	}

// Go over all couples
	std::size_t l = 0;
	double sum = 0;
	for (std::size_t i = 0; i < nbNodes; i++) {
		for (std::size_t j = i + 1; j < nbNodes; j++, l++) {
			double ksi = sqDist[l]; // Observed distance
			double w = weight[l];

			if (w > 0) { //i.e. known couple
				double d = 0.0;
				for (std::size_t k = 0; k < dim; k++) { // Compute the current distance between the two nodes
					std::size_t indexI = i * dim + k;
					std::size_t indexJ = j * dim + k;
					double v = x[indexI] - x[indexJ];
					d += v * v;
				}
				double val = std::log(d / ksi);
				sum += w * val * val;
				for (std::size_t k = 0; k < dim; k++) {
					std::size_t indexI = i * dim + k;
					std::size_t indexJ = j * dim + k;

					double v = w * (2.0 * val) * (x[indexI] - x[indexJ])
							/ (d * nbKnownCouples);
					grad_f[indexI] += v;
					grad_f[indexJ] -= v;
				}
			}
		}
	}
	return (sum / nbKnownCouples);
}

void LogMDSCG::finalize(double * x, CG::INT n) {
//#pragma omp parallel for
	for (std::size_t i = 0; i < nbNodes * (nbNodes - 1) / 2; i++) {
		sqDist[i] *= scl;
	}

	double sqScl = std::sqrt(scl);
//#pragma omp parallel for
	for (std::size_t i = 0; i < nbNodes * dim; i++) {
		x[i] *= sqScl;
	}
}

} /* namespace LinkPred */
