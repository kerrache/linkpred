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

int main(int argc, char *argv[]) {

	std::vector<Vec> trInput;
	trInput.push_back( { 0.912145, 0.709983, 0.226475 });
	trInput.push_back( { 0.934958, 0.123857, 0.802411 });
	trInput.push_back( { 0.039990, 0.781305, 0.560989 });
	trInput.push_back( { 0.322438, 0.241671, 0.637029 });
	trInput.push_back( { 0.895175, 0.726442, 0.406118 });
	trInput.push_back( { 0.140349, 0.068158, 0.488275 });
	trInput.push_back( { 0.474313, 0.968052, 0.370530 });
	trInput.push_back( { 0.437717, 0.953002, 0.371601 });
	trInput.push_back( { 0.655664, 0.527321, 0.712499 });
	trInput.push_back( { 0.123821, 0.552098, 0.846477 });

	std::vector<bool> trOutput = { 0, 1, 0, 1, 0, 1, 0, 0, 1, 0 };

	std::vector<Vec> tsInput;
	tsInput.push_back( { 0.85568, 0.36109, 0.86532 });
	tsInput.push_back( { 0.13094, 0.61792, 0.80714 });
	tsInput.push_back( { 0.61693, 0.47719, 0.67608 });
	tsInput.push_back( { 0.47321, 0.57101, 0.10932 });
	tsInput.push_back( { 0.73278, 0.19042, 0.70569 });

	std::vector<bool> tsOutput = { 1, 0, 1, 0, 1 };

	std::vector<double> pred(5);

	auto classifier = std::make_shared<LogisticRegresser<>>(0.001, 777);

	classifier->learn(trInput.begin(), trInput.end(), trOutput.begin(),
			trOutput.end());
	classifier->predict(tsInput.begin(), tsInput.end(), pred.begin());

	std::cout << "Pred\tActual" << std::endl;
	for (int j = 0; j < 5; j++) {
		std::cout << pred[j] << "\t" << tsOutput[j] << std::endl;
	}

	return 0;
}
