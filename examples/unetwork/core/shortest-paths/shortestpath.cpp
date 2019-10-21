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
	if (argc != 2) {
		std::cerr << "Bad arguments\nUsage: " << argv[0] << " netFileName\n";
		exit(1);
	}
	auto net = UNetwork<>::read(std::string(argv[1]));
	auto length = net->template createEdgeMapSP<double>();
	int i = 1;
	for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it, i++) {
		(*length)[*it] = (13 * i) % 3 + 1;
	}
	Dijkstra<> dijkstra(net);
	auto lengthMapId = dijkstra.registerLengthMap(length);
	{
		auto res = dijkstra.getShortestPath(0, net->getNbNodes() - 1,
				lengthMapId);
		auto path = res.first;
		auto dist = res.second;
		std::cout << "Path from " << net->getLabel(0) << " to "
				<< net->getLabel(net->getNbNodes() - 1) << " : ";
		for (auto it = path->begin(); it != path->end(); ++it) {
			std::cout << net->getLabel(*it) << ", ";
		}
		std::cout << "\nTotal distance: " << dist << std::endl;
	}
	{
		auto distMap = dijkstra.getDist(0, lengthMapId);
		std::cout << "Distance from node " << net->getLabel(0)
				<< " to all nodes:" << std::endl;
		std::cout << "Node\tDist\t# of hops\n";
		for (auto it = net->nodesBegin(); it != net->nodesEnd(); ++it) {
			auto res = distMap->at(it->first);
			std::cout << it->second << "\t" << res.first << "\t" << res.second
					<< std::endl;
		}
	}
	return 0;
}
