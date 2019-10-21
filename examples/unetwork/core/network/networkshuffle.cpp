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

int main(int argc, char*argv[]) {
	if (argc != 4) {
		std::cerr << "Bad arguments\nUsage:" << argv[0]
				<< " netFileName nbTests seed" << std::endl;
		exit(1);
	}
	std::string netFileName(argv[1]);
	int nbTests = std::atoi(argv[2]);
	long int seed = std::atol(argv[3]);

//	RandomGen rng(seed);
//	int n = 12;
//	std::vector<int> freq;
//	freq.resize(n, 0);
//	for (int i = 0; i < nbTests; i++) {
//		auto perm = Utilities::getRndPerm(n, rng.getInt());
//		for (std::size_t k = 0; k < n; k++) {
//			freq[k] += perm[k];
//		}
//	}
//	for (std::size_t k = 0; k < n; k++) {
//		std::cout << (double) freq[k] / nbTests << std::endl;
//	}

	auto net = UNetwork<>::read(netFileName, false, true);
	RandomGen rng(seed);

//	std::map<int, int> freq;
//	for (auto it = net->nodesBegin(); it != net->nodesEnd(); ++it) {
//		freq[it->second] = 0;
//	}
//
//	auto start = std::chrono::steady_clock::now();
//
//	for (int i = 0; i < nbTests; i++) {
//		net->shuffle(rng.getInt());
//		int k = 0;
//		for (auto it = net->nodesBegin(); it != net->nodesEnd(); ++it, k++) {
//			freq[it->second] += k;
//		}
//	}
//
//	for (auto it = net->nodesBegin(); it != net->nodesEnd(); ++it) {
//		std::cout << (double) freq[it->second] / nbTests << std::endl;
//	}

//	std::map<std::pair<int, int>, int> freq;
//	for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it) {
//		auto i = net->getLabel(net->start(*it));
//		auto j = net->getLabel(net->end(*it));
//		auto e = std::make_pair(std::min(i, j), std::max(i, j));
//		freq[e] = 0;
//	}
//
//	auto start = std::chrono::steady_clock::now();
//
//	for (int i = 0; i < nbTests; i++) {
//		net->shuffle(rng.getInt());
//		int k = 0;
//		for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it, k++) {
//			auto i = net->getLabel(net->start(*it));
//			auto j = net->getLabel(net->end(*it));
//			auto e = std::make_pair(std::min(i, j), std::max(i, j));
//			freq[e] = freq.at(e) + k;
//		}
//	}
//
//	for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it) {
//		auto i = net->getLabel(net->start(*it));
//		auto j = net->getLabel(net->end(*it));
//		auto e = std::make_pair(std::min(i, j), std::max(i, j));
//
//		std::cout << (double) freq.at(e) / nbTests << std::endl;
//	}
//
//	auto end = std::chrono::steady_clock::now();
//	auto diff = end - start;
//	std::cerr << "#Time: "
//			<< std::chrono::duration<double, std::milli>(diff).count() << " ms"
//			<< std::endl;

	std::map<std::pair<typename UNetwork<>::LabelType, typename UNetwork<>::LabelType>, int> freq;
	for (auto it = net->nonEdgesBegin(); it != net->nonEdgesEnd(); ++it) {
		auto i = net->getLabel(net->start(*it));
		auto j = net->getLabel(net->end(*it));
		auto e = std::make_pair(std::min(i, j), std::max(i, j));
		freq[e] = 0;
	}

	auto start = std::chrono::steady_clock::now();

	for (int i = 0; i < nbTests; i++) {
		net->shuffle(rng.getInt());
		int k = 0;
		for (auto it = net->nonEdgesBegin(); it != net->nonEdgesEnd();
				++it, k++) {
			auto i = net->getLabel(net->start(*it));
			auto j = net->getLabel(net->end(*it));
			auto e = std::make_pair(std::min(i, j), std::max(i, j));
			freq[e] = freq.at(e) + k;
		}
	}

	for (auto it = net->nonEdgesBegin(); it != net->nonEdgesEnd(); ++it) {
		auto i = net->getLabel(net->start(*it));
		auto j = net->getLabel(net->end(*it));
		auto e = std::make_pair(std::min(i, j), std::max(i, j));

		std::cout << (double) freq.at(e) / nbTests << std::endl;
	}

	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;
	std::cerr << "#Time: "
			<< std::chrono::duration<double, std::milli>(diff).count() << " ms"
			<< std::endl;

	return 0;
}
