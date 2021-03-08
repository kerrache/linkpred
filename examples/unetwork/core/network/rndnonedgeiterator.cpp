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
#include <map>

using namespace LinkPred;

int main(int argc, char*argv[]) {
	if (argc != 5) {
		std::cerr << "Bad arguments\nUsage:" << argv[0]
				<< " netFileName ratio nbTests seed" << std::endl;
		exit(1);
	}
	std::string netFileName(argv[1]);
	double ratio = std::atof(argv[2]);
	int nbTests = std::atoi(argv[3]);
	long int seed = std::atol(argv[4]);

	if (ratio < 0 || ratio > 1) {
		throw std::invalid_argument(
				"Illegal value for ratio: ratio must satisfy 0 <= ratio <= 1");
	}
	std::cout << "Reading network...\n";
	auto net = UNetwork<>::read(netFileName, false, true);
	std::cout << "done.\n";
	std::cout << "n: " << net->getNbNodes() << " m: " << net->getNbEdges()
			<< std::endl;
	net->print();

	{
		RandomGen rng(seed);

		std::map<typename UNetwork<>::Edge, int> freq;
		for (auto it = net->nonEdgesBegin(); it != net->nonEdgesEnd(); ++it) {
			freq[*it] = 0;
		}

		for (int i = 0; i < nbTests; i++) {
			for (auto it = net->rndNonEdgesBegin(ratio, rng.getInt());
					it != net->rndNonEdgesEnd(); ++it) { // pre-increment
				auto e = *it;
				auto fit = freq.find(e);
				if (fit == freq.end()) {
					throw std::logic_error("Unknown negative link");
				} else {
					fit->second++;
				}
			}
		}
		for (auto it = freq.begin(); it != freq.end(); ++it) {
			std::cout << net->getLabel(net->start(it->first)) << "\t"
					<< net->getLabel(net->end(it->first)) << "\t"
					<< (double) it->second / nbTests << std::endl;
		}
	}

	{
		RandomGen rng(seed);

		std::map<typename UNetwork<>::Edge, int> freq;
		for (auto it = net->nonEdgesBegin(); it != net->nonEdgesEnd(); ++it) {
			freq[*it] = 0;
		}

		for (int i = 0; i < nbTests; i++) {
			for (auto it = net->rndNonEdgesBegin(ratio, rng.getInt());
					it != net->rndNonEdgesEnd(); it++) { // post-increment
				auto e = *it;
				auto fit = freq.find(e);
				if (fit == freq.end()) {
					throw std::logic_error("Unknown negative link");
				} else {
					fit->second++;
				}
			}
		}
		for (auto it = freq.begin(); it != freq.end(); ++it) {
			std::cout << net->getLabel(net->start(it->first)) << "\t"
					<< net->getLabel(net->end(it->first)) << "\t"
					<< (double) it->second / nbTests << std::endl;
		}
	}
	return 0;
}
