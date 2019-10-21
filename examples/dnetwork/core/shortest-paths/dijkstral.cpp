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
int main(int argc, char*argv[]) {
	if (argc != 3) {
		std::cerr << "Bad arguments\nUsage: " << argv[0] << " netFileName lim\n";
		exit(1);
	}
	auto net = DNetwork<>::read(std::string(argv[1]), false, true);
	std::size_t lim = std::atol(argv[2]);
	auto length = net->template createEdgeMapSP<double>();
	int i = 1;
	std::cout << "Src\tEnd\tLength" << std::endl;
	for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it, i++) {
		//double l = (13 * i) % 3 + 1;
		double l = net->getOutDeg(net->start(*it)) + net->getInDeg(net->end(*it));
		(*length)[*it] = l;
		std::cout << net->getLabel(net->start(*it)) << "\t"
				<< net->getLabel(net->end(*it)) << "\t" << l << std::endl;
	}
	std::cout << "-------------------------------------\n";
	Dijkstra<DNetwork<>> dijkstra(net);
	auto lengthMapId = dijkstra.registerLengthMap(length);
	std::cout << "#Src\tDst\tDist\tnbHops" << std::endl;
	for (auto sit = net->nodesBegin(); sit != net->nodesEnd(); ++sit) {
		auto dist = dijkstra.getDistL(sit->first, lengthMapId, lim);
		for (auto dit = net->nodesBegin(); dit != net->nodesEnd(); ++dit) {
			auto res = dist->at(dit->first);
			std::cout << sit->second << "\t" << dit->second << "\t"
					<< res.first << "\t" << res.second << std::endl;
		}
	}
	return 0;
}
