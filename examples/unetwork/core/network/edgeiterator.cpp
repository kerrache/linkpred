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
	UNetwork<unsigned int> net;
	for (int i = 1; i <= n; i++) {
		net.addEdge(net.addNode(i).first, net.addNode(i % n + 1).first);
	}
	net.assemble();

	std::cout << "Positive links:" << std::endl;
	std::cout << "Start\tEnd" << std::endl;
	for (auto it = net.nodesDegBegin(); it != net.nodesDegEnd(); ++it) {
		for (auto nit = net.neighbBegin(it->first);
				nit != net.neighbEnd(it->first); ++nit) {
			std::cout << net.getLabel(net.start(*nit)) << "\t"
					<< net.getLabel(net.end(*nit)) << std::endl;
		}
	}

	std::cout << "Random negative links:" << std::endl;
	double ratio = 0.2;
	long int seed = 777;
	std::cout << "Start\tEnd" << std::endl;
	for (auto it = net.rndNonEdgesBegin(ratio, seed);
			it != net.rndNonEdgesEnd(); ++it) {
		std::cout << net.getLabel(net.start(*it)) << "\t"
				<< net.getLabel(net.end(*it)) << std::endl;
	}
	return 0;
}
