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
	std::cout << 0 << "\t";
	for (auto it = fullNet->nodesDegBegin(); it != fullNet->nodesDegEnd();
			++it) {
		std::cout << it->second << "\t";
	}
	std::cout << std::endl;

	for (std::size_t i = 0; i < nbTests; i++) {
		for (double ratio = 0.1; ratio <= 0.9; ratio += 0.1) {
			auto testData = NetworkManipulator<>::createTestData(fullNet, 0,
					ratio * fullNet->getNbEdges()
							/ (fullNet->getNbNodes()
									* (fullNet->getNbNodes() - 1) / 2), false,
					true, 0, true, 0, TP, FP, rng.getInt());
			std::cout << ratio << "\t";
			for (auto it = testData.getObsNet()->nodesDegBegin();
					it != testData.getObsNet()->nodesDegEnd(); ++it) {
				std::cout << it->second << "\t";
			}
			std::cout << std::endl;
		}
	}

	return 0;
}
