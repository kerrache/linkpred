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

#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;

int main(int argc, char *argv[]) {
	int nbNodes = 10;
	int dim = 2;
	double lambda = 0.0001;
	double tol = 1.0e-8;
	long int seed = 777;
	int n = nbNodes * dim;

	std::vector<double> trueCoords(n);
	RandomGen rng(seed);
	for (int i = 0; i < n; i++) {
		trueCoords[i] = rng.getDouble(0, 1);
	}
	std::vector<MatFactPbData> pbData;
	for (int i = 0; i < nbNodes; i++) {
		for (int j = i + 1; j < nbNodes; j++) {
			double target = 0;
			for (int k = 0; k < dim; k++) {
				target += trueCoords[i * dim + k] * trueCoords[j * dim + k];
			}
			pbData.push_back( { i, j, target });
		}
	}

	std::vector<double> coords(n);
	CG::cg_stats stats;
	CG::CGDProblem* pb = new MatFactCG(nbNodes, pbData, dim, lambda,
			rng.getInt());
	CG::CGDescent solver(pb);
	solver.cg_descent(coords.data(), n, &stats, nullptr, tol, nullptr, true);

	double obj = stats.f;
	std::cout << "obj: " << obj << std::endl;
	double rmse = 0;
	int k = 0;
	for (int i = 0; i < nbNodes; i++) {
		for (int j = i + 1; j < nbNodes; j++) {
			double actual = 0;
			for (int k = 0; k < dim; k++) {
				actual += coords[i * dim + k] * coords[j * dim + k];
			}
			double err = actual - pbData[k++].target;
			rmse += err * err;
		}
	}
	rmse = std::sqrt(rmse / pbData.size());
	std::cout << "rmse: " << rmse << std::endl;
	return 0;
}
