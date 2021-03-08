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
	std::cout << "Loading test data..." << std::endl;
	auto testData = NetworkManipulator<>::loadTestData("net-obs.edges",
			"net-rem.edges", "", true, 0, true, 0, FN, TN, 777, true);
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
