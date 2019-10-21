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
	if (argc != 2) {
		std::cerr << "Bad arguments\nUsage: " << argv[0] << " net" << std::endl;
		exit(1);
	}
	std::string netFileName(argv[1]);
	std::cout << "Reading network..." << std::endl;
	auto net = UNetwork<>::read(netFileName, false, true);
	auto start = std::chrono::steady_clock::now();
	std::cout << "Creating test set..." << std::endl;
	auto testData = NetworkManipulator<>::createTestData(net, 0.1, 0.0, false,
			true, 0, true, 0, FN, TN, 777, false);
	testData.genPos();
	testData.genNeg();
	std::cout << "Positive links in the test set (false negatives):"
			<< std::endl;
	int cpt = 0;
	for (auto it = testData.posBegin(); it != testData.posEnd(); ++it) {
		//std::cout << net->getLabel(net->start(*it)) << "\t" << net->getLabel(net->end(*it)) << std::endl;
		cpt++;
	}
	std::cout << "# pos edges: " << cpt << std::endl;
	std::cout << "Negative links in the test set (true negatives):"
			<< std::endl;
	cpt = 0;
	for (auto it = testData.negBegin(); it != testData.negEnd(); ++it) {
		//std::cout << net->getLabel(net->start(*it)) << "\t" << net->getLabel(net->end(*it)) << std::endl;
		cpt++;
	}
	std::cout << "# neg edges: " << cpt << std::endl;
	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;
	std::cerr << "#Time: "
			<< std::chrono::duration<double, std::milli>(diff).count() << " ms"
			<< std::endl;

	return 0;
}
