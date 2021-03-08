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
	int n = 8;
	DNetwork<unsigned int> net;
	for (int i = 1; i <= n; i++) {
		net.addEdge(net.addNode(i).first, net.addNode(i % n + 1).first);
		net.addEdge(net.addNode(i).first, net.addNode((i + 1) % n + 1).first);
	}
	net.assemble();

	int i = 0;
	auto nodeMap = net.template createNodeMap<double>();
	for (auto it = net.nodesBegin(); it != net.nodesEnd(); ++it) {
		nodeMap[it->first] = i++ / 2.0;
	}

	std::cout << "ID\tValue" << std::endl;
	for (auto it = net.nodesBegin(); it != net.nodesEnd(); ++it) {
		std::cout << it->second << "\t" << nodeMap.at(it->first) << std::endl;
	}

	auto edgeMap = net.template createEdgeMap<double>();
	i = 0;
	for (auto it = net.edgesBegin(); it != net.edgesEnd(); ++it) {
		double w = i++ / 2.0;
		edgeMap[*it] = w;
	}

	std::cout << "Start\tEnd\tValue" << std::endl;
	for (auto it = net.edgesBegin(); it != net.edgesEnd(); ++it) {
		std::cout << net.getLabel(net.start(*it)) << "\t"
				<< net.getLabel(net.end(*it)) << "\t" << edgeMap.at(*it)
				<< std::endl;
	}
	return 0;
}
