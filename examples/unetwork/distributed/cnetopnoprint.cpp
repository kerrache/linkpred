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

#include "linkpred.hpp"
#include <iostream>
#include <algorithm>
#include <chrono>

using namespace LinkPred;

int main(int argc, char*argv[]) {
#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
#endif
	if (argc != 3) {
		std::cerr << "Bad arguments\nUsage: " << argv[0]
				<< " netFileName top\n";
		exit(1);
	}
	auto start = std::chrono::steady_clock::now();
	std::string netFileName(argv[1]);
	std::size_t k = std::atol(argv[2]);

	auto net = UNetwork<>::read(netFileName, false, true);
	UCNEPredictor<> predictor(net);
#ifdef WITH_MPI
	predictor.setComm(MPI_COMM_WORLD); // Optional when using the default communicator MPI_COMM_WORLD
	predictor.setDistributed(true);
#endif
#ifdef WITH_OPENMP
	predictor.setParallel(true);
#endif

	predictor.init();
	predictor.learn();

	std::vector<typename UNetwork<>::EdgeType> edges;
	edges.resize(k);
	std::vector<double> scores;
	scores.resize(k);
	k = predictor.top(k, edges.begin(), scores.begin());

	int procID = 0;
#ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);
#endif

	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;
	if (procID == 0) {
		std::cerr << "#Time: "
				<< std::chrono::duration<double, std::milli>(diff).count()
				<< " ms" << std::endl;
	}
#ifdef WITH_MPI
	MPI_Finalize();
#endif
	return 0;
}
