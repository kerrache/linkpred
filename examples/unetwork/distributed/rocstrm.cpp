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
#include <chrono>
using namespace LinkPred;
int main(int argc, char *argv[]) {
#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
#endif
	if (argc != 2) {
		std::cerr << "Bad arguments\nUsage: " << argv[0] << " netFileName\n";
		exit(1);
	}
	int procID = 0;
#ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);
#endif

	std::string netFileName(argv[1]);
	auto fullNet = UNetwork<>::read(netFileName, false, true);
	auto testData = NetworkManipulator<>::createTestData(fullNet, 0.1, 0, false,
			true, 0, true, 0, FN, TN, 777, false);
	testData.lock();
	auto predictor = std::make_shared<UADAPredictor<>>(testData.getObsNet());
	predictor->init();
	predictor->learn();
	auto predResults = std::make_shared<PredResults<>>(testData, predictor);
	auto roc = std::make_shared<ROC<>>("ROC");
#ifdef WITH_OPENMP
	roc->setParallel(true);
#endif
#ifdef WITH_MPI
	roc->setComm(MPI_COMM_WORLD); // Optional when using the default communicator MPI_COMM_WORLD
	roc->setDistributed(true);
#endif
	roc->setStrmEnabled(true);
	PerfResults res;
	auto start = std::chrono::steady_clock::now();
	roc->eval(predResults, res);
	if (procID == 0) {
		std::cout << "#ROCAUC (streaming): " << res.at(roc->getName())
				<< std::endl;
		auto end = std::chrono::steady_clock::now();
		auto diff = end - start;
		std::cerr << "#Time: "
				<< std::chrono::duration<double, std::milli>(diff).count()
				<< " ms" << std::endl;
	}
#ifdef WITH_MPI
	MPI_Finalize();
#endif
	return 0;
}
