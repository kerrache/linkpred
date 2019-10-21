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
	if (argc != 4) {
		std::cerr << "Bad arguments\nUsage: " << argv[0]
				<< " netFileName top check\n";
		exit(1);
	}
	auto start = std::chrono::steady_clock::now();
	std::string netFileName(argv[1]);
	std::size_t k = std::atol(argv[2]);
	bool check = (std::atoi(argv[3]) != 0);

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

	std::cout << "Process " << procID << "\t" << k << std::endl;
	if (procID == 0) {
		std::cout << "#Start\tEnd\tScore\n";
	}
	for (std::size_t i = 0; i < k; i++) {
		std::cout << net->getLabel(net->start(edges[i])) << "\t"
				<< net->getLabel(net->end(edges[i])) << "\t" << scores[i]
				<< std::endl;
	}

	// Checking results if requested
	if (procID == 0 && check) {
		std::vector<typename UNetwork<>::EdgeType> edges;
		edges.resize(net->getNbNonEdges());
		std::copy(net->nonEdgesBegin(), net->nonEdgesEnd(), edges.begin());
		std::vector<double> scores;
		scores.resize(net->getNbNonEdges());
		predictor.predict(edges.begin(), edges.end(), scores.begin());

		LMapQueue<typename UNetwork<>::EdgeType, double> mq(k);
		std::size_t i = 0;
		for (auto it = edges.begin(); it != edges.end(); ++it, i++) {
			mq.push(*it, scores[i]);
		}

		std::cout << "#Start\tEnd\tScore\n";
		i = 0;
		for (auto it = mq.begin(); it != mq.end(); ++it, i++) {
			std::cout << net->getLabel(net->start(it->first)) << "\t"
					<< net->getLabel(net->end(it->first)) << "\t" << it->second
					<< std::endl;
		}
	}

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
