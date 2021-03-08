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
#include <algorithm>
using namespace LinkPred;

int main(int argc, char*argv[]) {
	if (argc != 2) {
		std::cerr << "Bad arguments\nUsage: " << argv[0] << " netFileName\n";
		exit(1);
	}
#ifdef LINKPRED_WITH_OPENMP
	omp_set_nested(1); // enable nested parallelism
#endif
	std::string netFileName(argv[1]);
	auto net = UNetwork<>::read(netFileName, false, true);
	UCRAPredictor<> predictor(net);
#ifdef LINKPRED_WITH_OPENMP
	predictor.setParallel(true);
#endif
	predictor.init();
	predictor.learn();

	std::vector<typename UNetwork<>::Edge> edges;
	edges.resize(net->getNbNonEdges());
	std::copy(net->nonEdgesBegin(), net->nonEdgesEnd(), edges.begin());
	std::vector<double> scores;
	scores.resize(net->getNbNonEdges());
	predictor.predict(edges.begin(), edges.end(), scores.begin());

	/*
	 std::cout << "#Start\tEnd\tScore\n";
	 std::size_t i = 0;
	 for (auto it = edges.begin(); it != edges.end(); ++it, i++) {
	 std::cout << net->getLabel(net->start(*it)) << "\t"
	 << net->getLabel(net->end(*it)) << "\t" << scores[i]
	 << std::endl;
	 }
	 */
	return 0;
}
