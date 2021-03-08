/*
 * This file is part of LinkPred.
 *
 * LinkPred: A high performance library for link prediction in complex networks.
 * Copyright (C) 2021  by Said Kerrache.
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
#ifdef LINKPRED_WITH_MLPACK

	arma::mat trainInput;
	mlpack::data::Load("train-input.csv", trainInput, true);
	arma::mat trainOutput;
	mlpack::data::Load("train-output.csv", trainOutput, true);
	arma::mat testInput;
	mlpack::data::Load("test-input.csv", testInput, true);
	arma::mat testOutput;
	mlpack::data::Load("test-output.csv", testOutput, true);

	int m = trainInput.n_rows; // Number of features

	FFN<> classifier;
	classifier.setAutoArch(m); // Automatic architecture

	// Marshalling
	int n = trainInput.n_cols; // Number of examples in the training set
	std::vector<Vec> trInput;
	for (int j = 0; j < n; j++) {
		Vec v(m);
		for (int i = 0; i < m; i++) {
			v[i] = trainInput(i, j);
		}
		trInput.push_back(v);
	}

	std::vector<bool> trOutput;
	for (int j = 0; j < n; j++) {
		trOutput.push_back(trainOutput(j));
	}

	// Training
	classifier.learn(trInput.begin(), trInput.end(), trOutput.begin(),
			trOutput.end());

	// Marshalling
	n = testInput.n_cols; // Number of examples in the test set
	std::vector<Vec> tsInput;
	for (int j = 0; j < n; j++) {
		Vec v(m);
		for (int i = 0; i < m; i++) {
			v[i] = testInput(i, j);
		}
		tsInput.push_back(v);
	}

	std::vector<double> pred;
	pred.resize(n);

	// Predicting
	classifier.predict(tsInput.begin(), tsInput.end(), pred.begin());

	// Computing error
	double rmse = 0;
	double mae = 0;
	for (int j = 0; j < n; j++) {
		double err = std::abs(pred[j] - testOutput[j]);
		mae += err;
		rmse += err * err;
	}
	rmse = std::sqrt(rmse / n);
	mae = mae / n;
	std::cout << "rmse: " << rmse << std::endl;
	std::cout << "mae : " << mae << std::endl;

#endif
	return 0;
}
