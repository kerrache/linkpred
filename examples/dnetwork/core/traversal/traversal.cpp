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
		std::cerr << "Bad arguments\nUsage: " << argv[0] << " netFileName srcNode\n";
		exit(1);
	}
	auto net = DNetwork<>::read(std::string(argv[1]));
	std::string srcNode = argv[2];
	// BFS
	BFS<DNetwork<>> bfs(net);
	Collector<DNetwork<>> col;
	bfs.traverse(net->getID(srcNode), col);
	auto visited = col.getVisited();
	std::cout << "BFS:" << std::endl;
	while (!visited.empty()) {
		auto i = visited.front();
		visited.pop();
		std::cout << net->getLabel(i) << std::endl;
	}
	// DFS
	DFS<DNetwork<>, Counter<DNetwork<>>> dfs(net);
	Counter<DNetwork<>> counter;
	dfs.traverse(net->getID(srcNode), counter);
	std::cout << "DFS visited " << counter.getCount() << " nodes" << std::endl;
	return 0;
}

