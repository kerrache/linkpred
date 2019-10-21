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
	bool first = true;
	for (std::size_t i = 0; i < nbTests; i++) {
		for (double ratio = 0.1; ratio <= 0.9; ratio += 0.1) {
			auto testData = NetworkManipulator<>::createTestData(fullNet, ratio,
					0, false, true, 0, true, 0, FN, TN, rng.getInt());
			testData.lock();
			PerfEvaluator<> perf(testData);

			{
				auto predictor = std::make_shared<UADAPredictor<>>(
						testData.getObsNet());
#ifdef WITH_OPENMP
				predictor->setParallel(true);
#endif
				perf.addPredictor(predictor);
			}
			{
				auto predictor = std::make_shared<UCNEPredictor<>>(
						testData.getObsNet());
#ifdef WITH_OPENMP
				predictor->setParallel(true);
#endif
				perf.addPredictor(predictor);
			}
			{
				auto predictor = std::make_shared<USHPPredictor<>>(
						testData.getObsNet(), rng.getInt());
#ifdef WITH_OPENMP
				predictor->setParallel(true);
#endif
				perf.addPredictor(predictor);
			}
			perf.addPerfMeasure(
					std::make_shared<TPR<>>(
							fullNet->getNbEdges()
									- testData.getObsNet()->getNbEdges()));
			perf.addPerfMeasure(std::make_shared<PR<>>());
			perf.addPerfMeasure(std::make_shared<ROC<>>());
			perf.eval();

			if (first) {
				std::cout << "#ratio\t";
				for (auto it = perf.resultsBegin(); it != perf.resultsEnd();
						++it) {
					std::cout << it->first << "\t";
				}
				std::cout << std::endl;
				first = false;
			}
			std::cout << std::fixed << std::setprecision(2) << ratio << "\t";
			for (auto it = perf.resultsBegin(); it != perf.resultsEnd(); ++it) {
				std::cout << std::fixed << std::setprecision(4) << it->second
						<< "\t";
			}
			std::cout << std::endl;
		}
	}
	return 0;
}
