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

int main(int argc, char *argv[]) {
	int n = 6;
	auto net = std::make_shared<UNetwork<unsigned int>>();
	for (int i = 1; i <= n; i++) {
		net->addEdge(net->addNode(i).first, net->addNode(i % n + 1).first);
		net->addEdge(net->addNode(i).first,
				net->addNode((i + 1) % n + 1).first);
	}
	net->assemble();
	std::cout << "Assembled\n";
	auto testData = NetworkManipulator<UNetwork<unsigned int>>::createTestData(
			net, 0.4, 0.3, true, true, 0, true, 0, FN, TN, 777);
	testData.genPos();
	testData.genNeg();
	std::cout << "Reference (original) network:\n";
	testData.getRefNet()->print();
	std::cout << "Observed (modified) network:\n";
	testData.getObsNet()->print();
	std::cout << "Positive links in the test set (false negatives):"
			<< std::endl;
	for (auto it = testData.posBegin(); it != testData.posEnd(); ++it) {
		std::cout << net->getLabel(net->start(*it)) << "\t"
				<< net->getLabel(net->end(*it)) << std::endl;
	}
	std::cout << "Negative links in the test set (true negatives):"
			<< std::endl;
	for (auto it = testData.negBegin(); it != testData.negEnd(); ++it) {
		std::cout << net->getLabel(net->start(*it)) << "\t"
				<< net->getLabel(net->end(*it)) << std::endl;
	}
	return 0;
}
