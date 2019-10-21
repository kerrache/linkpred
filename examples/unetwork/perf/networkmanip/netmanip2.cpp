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
int main(int argc, char *argv[]) {
	if (argc != 3) {
		std::cerr << "Bad arguments\nUsage: " << argv[0] << " net1 net2"
				<< std::endl;
		exit(1);
	}
	std::string net1FileName(argv[1]);
	std::string net2FileName(argv[2]);
	std::cout << "Reading networks..." << std::endl;
	auto net1 = UNetwork<>::read(net1FileName, false, true);
	auto net2 = UNetwork<>::read(net2FileName, false, true);
	std::cout << "Creating test set..." << std::endl;
	auto testData = NetworkManipulator<>::createTestDataSeq(net1, net2, true, 0,
			true, 0, FN, TN, 777);
	testData.genPos();
	testData.genNeg();
	std::cout << "First network:\n";
	net1->print();
	std::cout << "Second network:\n";
	net2->print();
	std::cout << "Test data original network:\n";
	testData.getRefNet()->print();
	std::cout << "Test data modified network:\n";
	testData.getObsNet()->print();
	std::cout << "Positive links in the test set (false negatives):"
			<< std::endl;
	for (auto it = testData.posBegin(); it != testData.posEnd(); ++it) {
		std::cout
				<< testData.getRefNet()->getLabel(
						testData.getRefNet()->start(*it)) << "\t"
				<< testData.getRefNet()->getLabel(
						testData.getRefNet()->end(*it)) << std::endl;
	}
	std::cout << "Negative links in the test set (true negatives):"
			<< std::endl;
	for (auto it = testData.negBegin(); it != testData.negEnd(); ++it) {
		std::cout
				<< testData.getRefNet()->getLabel(
						testData.getRefNet()->start(*it)) << "\t"
				<< testData.getRefNet()->getLabel(
						testData.getRefNet()->end(*it)) << std::endl;
	}
	return 0;
}
