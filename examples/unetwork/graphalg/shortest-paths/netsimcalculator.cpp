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
	auto length = net->template createEdgeMapSP<double>();
	int i = 1;
	std::cout << "Src\tEnd\tLength" << std::endl;
	for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it, i++) {
		//double l = (13 * i) % 3 + 1;
		double l = (net->getDeg(net->start(*it)) + net->getDeg(net->end(*it)))
				/ (2.0 * net->getMaxDeg());
		(*length)[*it] = l;
		std::cout << net->getLabel(net->start(*it)) << "\t"
				<< net->getLabel(net->end(*it)) << "\t" << l << std::endl;
	}
	Dijkstra<> dijkstra(net);
	{
		ESPSimlCalculator<> calc(dijkstra, length, NoCache);
		std::cout << "Src\tDst\tSiml\tDist\tnbHops" << std::endl;
		for (auto sit = net->nodesBegin(); sit != net->nodesEnd(); ++sit) {
			for (auto dit = net->nodesBegin(); dit != net->nodesEnd(); ++dit) {
				auto res = calc.getSiml(sit->first, dit->first);
				std::cout << sit->second << "\t" << dit->second << "\t"
						<< std::get<0>(res) << "\t" << std::get<1>(res) << "\t"
						<< std::get<2>(res) << "\t" << std::endl;
			}
		}
	}
	/*
	 {
	 ESPDsimCalculator<> calc(dijkstra, length, NoCache);
	 std::cout << "Src\tDst\tSiml\tDist\tnbHops" << std::endl;
	 for (auto sit = net->nodesBegin(); sit != net->nodesEnd(); ++sit) {
	 for (auto dit = net->nodesBegin(); dit != net->nodesEnd(); ++dit) {
	 auto resi = calc.getDist(sit->first, dit->first);
	 auto resj = calc.getDist(dit->first, sit->first);
	 std::cout << sit->second << "\t" << dit->second << "\t"
	 << 0.5 * (resi.first + resj.first) << "\t"
	 << resi.second << std::endl;
	 }
	 }
	 }
	 */
	return 0;
}
