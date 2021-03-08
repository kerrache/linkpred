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

int main(int argc, char*argv[]) {
	if (argc != 2) {
		std::cerr << "Bad arguments\nUsage: " << argv[0] << " netFileName\n";
		exit(1);
	}
	auto net = UNetwork<>::read(std::string(argv[1]));
	net->print();
	auto length = net->template createEdgeMapSP<double>();
	int i = 1;
	for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it, i++) {
		(*length)[*it] = (13 * i) % 3 + 1;
	}
	Dijkstra<> dijkstra(net);
	std::size_t lim = 2;
	ESPLDistCalculator<> calc(dijkstra, length, lim, NetworkCache);
	std::cout << "Src\tDst\tDist\tnbHops" << std::endl;
	for (auto sit = net->nodesBegin(); sit != net->nodesEnd(); ++sit) {
		for (auto dit = sit + 1; dit != net->nodesEnd(); ++dit) {
			auto dst = calc.getDist(sit->first, dit->first);
			std::cout << sit->second << "\t" << dit->second << "\t" << dst.first
					<< "\t" << dst.second << std::endl;
		}
	}
	return 0;
}
