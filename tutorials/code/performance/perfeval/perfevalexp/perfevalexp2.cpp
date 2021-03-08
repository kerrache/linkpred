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

using namespace LinkPred;

class Factory: public PEFactory<> {
public:
	virtual std::vector<std::shared_ptr<ULPredictor<>>> getPredictors(
			std::shared_ptr<UNetwork<> const> obsNet) {
		std::vector<std::shared_ptr<ULPredictor<>>> prs;
		prs.push_back(std::make_shared<URALPredictor<>>(obsNet));
		prs.push_back(std::make_shared<UKABPredictor<>>(obsNet));
		return prs;
	}
	virtual std::vector<std::shared_ptr<PerfMeasure<>>> getPerfMeasures(
			TestData<> const & testData) {
		std::vector<std::shared_ptr<PerfMeasure<>>> pms;
		pms.push_back(std::make_shared<TPR<>>(testData.getNbPos()));
		pms.push_back(std::make_shared<ROC<>>());
		return pms;
	}
	virtual ~Factory() = default;
};

int main(int argc, char*argv[]) {
	if (argc != 2) {
		std::cerr << "Bad arguments\nUsage: " << argv[0] << " netFileName\n";
		exit(1);
	}
	std::string netFileName(argv[1]);
	auto refNet = UNetwork<>::read(netFileName);
	PerfeEvalExpDescp<> ped;
	ped.refNet = refNet;
	ped.nbTestRuns = 100;
	ped.seed = 777;
	auto factory = std::make_shared<Factory>();
	PerfEvalExp<> exp(ped, factory);
	exp.run();
	return 0;
}

