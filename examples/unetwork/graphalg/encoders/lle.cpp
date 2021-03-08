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

template<typename Iter> void print(Iter it, int n, std::string name) {
	std::cout << name << ": ";
	for (int i = 0; i < n; i++) {
		std::cout << std::fixed << std::setprecision(6) << it[i] << " ";
	}
	std::cout << std::endl;
}

int main() {
#ifdef LINKPRED_WITH_ARMADILLO
	int  k = 3;
	int dim = 2;
	// Generate a regular k x k grid network
	auto net = UNetwork<>::generateREG(k, k);
	LLE<UNetwork<unsigned int>> encoder(net);
	encoder.setDim(dim);
	encoder.init();
	encoder.encode();
	for (std::size_t i = 0; i < net->getNbNodes(); i++) {
		auto v = encoder.getNodeCode(i);
		print(v, v.size(), "v");
	}
#endif
	return 0;
}
