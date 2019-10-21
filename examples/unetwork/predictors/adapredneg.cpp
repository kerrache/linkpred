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
#include <algorithm>
#include <chrono>

using namespace LinkPred;
int main(int argc, char*argv[]) {
	if (argc != 3) {
		std::cerr << "Bad arguments\nUsage: " << argv[0]
				<< " netFileName print[0/1]\n";
		exit(1);
	}
	auto start = std::chrono::steady_clock::now();
	std::string netFileName(argv[1]);
	bool print = (std::atoi(argv[2]) != 0);

	auto net = UNetwork<>::read(netFileName, false, true);
	UADAPredictor<> predictor(net);
	predictor.init();
	predictor.learn();

	std::vector<double> scores;
	scores.resize(net->getNbNonEdges());
	auto its = predictor.predictNeg(scores.begin());

	if (print) {
		std::cout << "#Start\tEnd\tScore\n";
		std::size_t i = 0;
		for (auto it = its.first; it != its.second; ++it, i++) {
			std::cout << net->getLabel(net->start(*it)) << "\t"
					<< net->getLabel(net->end(*it)) << "\t" << scores[i]
					<< std::endl;
		}
	}

	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;
	std::cerr << "#Time: "
			<< std::chrono::duration<double, std::milli>(diff).count()
			<< " ms" << std::endl;
	return 0;
}
