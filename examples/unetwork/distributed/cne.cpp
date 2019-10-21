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
using namespace LinkPred;

int main(int argc, char*argv[]) {
#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
#endif
	if (argc != 2) {
		std::cerr << "Bad arguments\nUsage: " << argv[0] << " netFileName\n";
		exit(1);
	}
	std::string netFileName(argv[1]);
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

	std::vector<double> scores;
#ifdef WITH_MPI
	scores.resize(predictor.localSize());
#else
	scores.resize(net->getNbNonEdges());
#endif
	auto range = predictor.predictNeg(scores.begin());

	int procID = 0;
#ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);
#endif
	std::cout << "#procID\tStart\tEnd\tScore\n";

	std::size_t i = 0;
	for (auto it = range.first; it != range.second; ++it, i++) {
		std::cout << procID << "\t";
		std::cout << net->getLabel(net->start(*it)) << "\t"
				<< net->getLabel(net->end(*it)) << "\t" << scores[i]
				<< std::endl;
	}

#ifdef WITH_MPI
	MPI_Finalize();
#endif
	return 0;
}
