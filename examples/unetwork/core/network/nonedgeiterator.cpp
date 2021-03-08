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
#include <chrono>

using namespace LinkPred;

int main(int argc, char*argv[]) {
	if (argc != 2) {
		std::cerr << "Bad arguments\nUsage:" << argv[0] << " netFileName"
				<< std::endl;
		exit(1);
	}
	std::string netFileName(argv[1]);
	std::cout << "Reading network...\n";
	auto net = UNetwork<>::read(netFileName, false, true);
	std::cout << "done.\n";
	std::cout << "n: " << net->getNbNodes() << " m: " << net->getNbEdges()
			<< std::endl;
	//net->print();
	//std::cout << "Negative links:\n";
	{
		auto start = std::chrono::steady_clock::now();
		std::size_t cpt = 0;
		for (auto it = net->nonEdgesBegin(); it != net->nonEdgesEnd(); ++it) { // pre-increment
			auto e = *it;
			if (e > 0) { // Just a dummy use of the variable e
				cpt++;
			}
			//std::cout << net->getLabel(net->start(e)) << "\t" << net->getLabel(net->end(e)) << std::endl;
		}
		auto end = std::chrono::steady_clock::now();
		auto diff = end - start;
		std::cerr << "#Time: "
				<< std::chrono::duration<double, std::milli>(diff).count()
				<< " ms" << std::endl;
		std::cout << cpt << " negative links " << std::endl;
	}
	//std::cout << "Negative links:\n";
	{
		auto start = std::chrono::steady_clock::now();
		std::size_t cpt = 0;
		for (auto it = net->nonEdgesBegin(); it != net->nonEdgesEnd(); it++) { // post-increment
			auto e = *it;
			if (e > 0) { // Just a dummy use of the variable e
				cpt++;
			}
			//std::cout << net->getLabel(net->start(e)) << "\t" << net->getLabel(net->end(e)) << std::endl;
		}
		auto end = std::chrono::steady_clock::now();
		auto diff = end - start;
		std::cerr << "#Time: "
				<< std::chrono::duration<double, std::milli>(diff).count()
				<< " ms" << std::endl;
		std::cout << cpt << " negative links " << std::endl;
	}
	{
		auto start = std::chrono::steady_clock::now();
		std::size_t cpt = 0;
		for (auto it = net->nonEdgesEnd(); it != net->nonEdgesBegin();) {
			it--;
			auto e = *it;
			if (e > 0) { // Just a dummy use of the variable e
				cpt++;
			}
			//std::cout << net->getLabel(net->start(e)) << "\t" << net->getLabel(net->end(e)) << std::endl;
		}
		auto end = std::chrono::steady_clock::now();
		auto diff = end - start;
		std::cerr << "#Time: "
				<< std::chrono::duration<double, std::milli>(diff).count()
				<< " ms" << std::endl;
		std::cout << cpt << " negative links " << std::endl;
	}
	{
		auto start = std::chrono::steady_clock::now();
		std::size_t cpt = 0;
		for (auto it = net->nonEdgesEnd(); it != net->nonEdgesBegin();) {
			--it;
			auto e = *it;
			if (e > 0) { // Just a dummy use of the variable e
				cpt++;
			}
			//std::cout << net->getLabel(net->start(e)) << "\t" << net->getLabel(net->end(e)) << std::endl;
		}
		auto end = std::chrono::steady_clock::now();
		auto diff = end - start;
		std::cerr << "#Time: "
				<< std::chrono::duration<double, std::milli>(diff).count()
				<< " ms" << std::endl;
		std::cout << cpt << " negative links " << std::endl;
	}
	return 0;
}
