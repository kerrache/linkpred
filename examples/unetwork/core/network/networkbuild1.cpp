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
	int n = 8;
	UNetwork<unsigned int> net;
	std::cout << "Label\tID\tNew?" << std::endl;
	for (int i = 1; i <= n; i++) {
		auto res = net.addNode(i);
		std::cout << i << "\t" << res.first << "\t" << res.second << std::endl;
	}
	for (int i = 1; i <= n; i++) {
		net.addEdge(net.getID(i), net.getID(i % n + 1));
		net.addEdge(net.getID(i), net.getID((i + 1) % n + 1));
	}
	net.assemble();
	std::cout << "Printing network:" << std::endl;
	net.print();
	return 0;
}
