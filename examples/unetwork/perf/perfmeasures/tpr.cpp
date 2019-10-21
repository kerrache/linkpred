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
	if (argc != 2) {
		std::cerr << "Bad arguments\nUsage: " << argv[0] << " netFileName\n";
		exit(1);
	}
	std::string netFileName(argv[1]);
	auto fullNet = UNetwork<>::read(netFileName, false, true);
	auto testData = NetworkManipulator<>::createTestData(fullNet, 0.3, 0, true,
			true, 0, true, 0, FN, TN, 777);
	auto predictor = std::make_shared<UHRGPredictor<>>(testData.getObsNet(),
			333);
	predictor->init();
	predictor->learn();
	auto predResults = std::make_shared<PredResults<>>(testData, predictor);
	auto tpr = std::make_shared<TPR<>>(100, "TPR");
	PerfResults res;
	tpr->eval(predResults, res);
	std::cout << "TPR: " << res.at(tpr->getName()) << std::endl;
	return 0;
}
