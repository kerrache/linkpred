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
	std::cout << "Reading networks..." << std::endl;
	auto net1 = UNetwork<>::read("net-seq1.edges");
	auto net2 = UNetwork<>::read("net-seq2.edges");
	std::cout << "First network:\n";
	net1->print();
	std::cout << "Second network:\n";
	net2->print();
	std::cout
			<< "Creating test set. Detecting links that will disappear: Positive class: TP. Negative class: FP"
			<< std::endl;
	auto testData = NetworkManipulator<>::createTestDataSeqInter(net1, net2,
			true, 0, true, 0, TP, FP, 777, true);
	std::cout << "Test data reference network:\n";
	auto refNet = testData.getRefNet();
	refNet->print();
	std::cout << "Test data observed network:\n";
	auto obsNet = testData.getObsNet();
	obsNet->print();
	std::cout << "Positive links in the test set:\n";
	for (auto it = testData.posBegin(); it != testData.posEnd(); ++it) {
		std::cout << refNet->getLabel(refNet->start(*it)) << "\t"
				<< refNet->getLabel(refNet->end(*it)) << std::endl;
	}
	std::cout << "Negative links in the test set:\n";
	for (auto it = testData.negBegin(); it != testData.negEnd(); ++it) {
		std::cout << refNet->getLabel(refNet->start(*it)) << "\t"
				<< refNet->getLabel(refNet->end(*it)) << std::endl;
	}
	return 0;
}
