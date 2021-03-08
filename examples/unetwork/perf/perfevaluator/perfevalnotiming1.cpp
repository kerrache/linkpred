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

#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;

int main(int argc, char*argv[]) {
	if (argc != 4) {
		std::cerr << "Bad arguments\nUsage: " << argv[0]
				<< " netFileName seed nbTests\n";
		exit(1);
	}

	std::string netFileName(argv[1]);
	long int seed = std::atol(argv[2]);
	std::size_t nbTests = std::atol(argv[3]);
	RandomGen rng(seed);
	auto fullNet = UNetwork<>::read(netFileName, false, true);
	for (std::size_t i = 0; i < nbTests; i++) {
		auto testData = NetworkManipulator<>::createTestData(fullNet, 0.1, 0,
				true, true, 0, true, 0, FN, TN, rng.getInt());
		PerfEvaluator<> perf(testData);

		perf.addPredictor(
				std::make_shared<UADAPredictor<>>(testData.getObsNet()));
		perf.addPredictor(
				std::make_shared<UCNEPredictor<>>(testData.getObsNet()));

		perf.addPerfMeasure(std::make_shared<ROC<>>());
		perf.addPerfMeasure(std::make_shared<PR<>>());

		perf.eval();

		if (i == 0) {
			std::cout << "#";
			for (auto it = perf.resultsBegin(); it != perf.resultsEnd(); ++it) {
				std::cout << it->first << "\t";
			}
			std::cout << std::endl;
		}
		for (auto it = perf.resultsBegin(); it != perf.resultsEnd(); ++it) {
			std::cout << std::setprecision(4) << it->second << "\t";
		}
		std::cout << std::endl;
	}
	return 0;
}
