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
		std::cerr << "Bad arguments\nUsage:" << argv[0] << " netFileName"
				<< std::endl;
		exit(1);
	}
	std::string netFileName(argv[1]);
	auto net = UNetwork<>::read(netFileName, false, true);

	netFileName = netFileName.substr(netFileName.find_last_of("/\\") + 1);
	std::cout << "#network\tn\tm\tminDeg\tmaxDeg\tavgDeg\tcc" << std::endl;
//	std::cout << net->getNbNodes() << "\t" << net->getNbEdges() << std::endl;
	std::cout << netFileName << "\t" << net->getNbNodes() << "\t"
			<< net->getNbEdges() << "\t" << net->getMinDeg() << "\t"
			<< net->getMaxDeg() << "\t" << net->getAvgDeg() << "\t"
			<< net->getCC() << std::endl;
	return 0;
}
