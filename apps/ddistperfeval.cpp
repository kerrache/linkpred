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
#include <fstream>
#include <vector>
#include <map>
#include <stdexcept>
#include <memory>
#include <chrono>
#include <sstream>
#include <iterator>

using namespace LinkPred;

/**
 * Method used to solve the link prediction problem.
 */
enum LPMethod {
	DCNE, /**< Directed common neighbors. */
	DADA,
	DSAI,
	DJID,
	DSOI,
	DHDI,
	DHPI,
	DPAT,
	DLHN,
	DLCP
};

std::map<std::string, LPMethod> lpNameMap = { { "DCNE", DCNE }, { "DADA", DADA }, { "DSAI", DSAI }, { "DJID", DJID }, { "DSOI", DSOI }, { "DHDI", DHDI }, { "DHPI", DHPI }, { "DPAT", DPAT }, { "DLHN", DLHN },{ "DLCP", DLCP } };

enum DSIPParams {
	Lim, UseSUM, Cache
};

std::map<std::string, DSIPParams> dsipParamsMap = { { "Lim", Lim }, { "UseSUM", UseSUM },
		{ "Cache", Cache } };

/**
 * Performance measures.
 */
enum PerfM {
	ROCT, PRT, TPRT, TPRTT
};

std::map<std::string, PerfM> pmNameMap = { { "ROC", ROCT }, { "PR", PRT }, {
		"TPR", TPRT }, { "TPRT", TPRTT } };

template<typename Out>
void split(const std::string &s, char delim, Out result) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		*(result++) = item;
	}
}

std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, std::back_inserter(elems));
	return elems;
}

std::shared_ptr<DLPredictor<>> createPredictor(std::string lpName,
		TestData<DNetwork<>> const & testData, RandomGen & rng) {

	auto tokens = split(lpName, '_');
	auto lpm = lpNameMap.at(tokens[0]);
	switch (lpm) {

	case DCNE:
		return std::make_shared<DCNEPredictor<>>(testData.getObsNet()); /**< Directed common neighbors. */

	case DADA:
		return std::make_shared<DADAPredictor<>>(testData.getObsNet());

	case DSAI:
		return std::make_shared<DSAIPredictor<>>(testData.getObsNet());

	case DJID:
		return std::make_shared<DJIDPredictor<>>(testData.getObsNet());

	case DSOI:
		return std::make_shared<DSOIPredictor<>>(testData.getObsNet());

	case DHDI:
		return std::make_shared<DHDIPredictor<>>(testData.getObsNet());

	case DHPI:
		return std::make_shared<DHPIPredictor<>>(testData.getObsNet());

	case DPAT:
		return std::make_shared<DPATPredictor<>>(testData.getObsNet());

	case DLHN:
		return std::make_shared<DLHNPredictor<>>(testData.getObsNet());

	case DLCP:
		return std::make_shared<DLCPPredictor<>>(testData.getObsNet());

	default:
		throw std::invalid_argument("Unknown link prediction method");
	}
}

void readFile(std::string fileName, bool & keepConnected, double & fnRatio,
		double &tnRatio, std::size_t & nbTests, double & ratioStart,
		double & ratioStep, double & ratioEnd,
		std::vector<PerfM> & perfMeasures,
		std::vector<std::string> & predictors, int procID) {

	std::string word;
	std::ifstream file(fileName);
	if (file.is_open()) {
		file >> keepConnected >> fnRatio >> tnRatio >> nbTests >> ratioStart
				>> ratioStep >> ratioEnd;
		if (procID == 0) {
			std::cout << "# keepConnected: " << keepConnected << std::endl;
			std::cout << "# fnRatio: " << fnRatio << std::endl;
			std::cout << "# tnRatio: " << tnRatio << std::endl;
			std::cout << "# nbTests: " << nbTests << std::endl;
			std::cout << "# ratioStart: " << ratioStart << std::endl;
			std::cout << "# ratioStep: " << ratioStep << std::endl;
			std::cout << "# ratioEnd: " << ratioEnd << std::endl;
		}
		while (file >> word) {
			auto fit = pmNameMap.find(word);
			if (fit != pmNameMap.end()) {
				perfMeasures.push_back(fit->second);
				if (procID == 0) {
					std::cout << "# Performance measure: " << word << std::endl;
				}
			} else {
				predictors.push_back(word);
				if (procID == 0) {
					std::cout << "# Predictor: " << word << std::endl;
				}
			}
		}
		file.close();
	} else {
		throw std::runtime_error("Unable to open file " + fileName);
	}
}

std::shared_ptr<PerfMeasure<PredResults<TestData<DNetwork<>>, DLPredictor<>>>> createPerfM(
		PerfM pm, long int nb) {
	switch (pm) {
	case ROCT:
		return std::make_shared<
				ROC<PredResults<TestData<DNetwork<>>, DLPredictor<>>>>();
	case PRT:
		return std::make_shared<
				PR<PredResults<TestData<DNetwork<>>, DLPredictor<>>>>();
	case TPRT:
		return std::make_shared<
				TPR<PredResults<TestData<DNetwork<>>, DLPredictor<>>>>(nb,
				"TPR");
	case TPRTT: {
		auto tpr = std::make_shared<
				TPR<PredResults<TestData<DNetwork<>>, DLPredictor<>>>>(nb,
				"TPRT");
		tpr->setUseTopMethod(true);
		return tpr;
	}
	default:
		throw std::invalid_argument("Unknown performance measure");
	}
}

int main(int argc, char*argv[]) {

	int procID = 0;
#ifdef WITH_MPI
	bool distributed = true;
	int nbProcs = 1;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nbProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);
#endif
#ifdef WITH_OPENMP
	bool parallel = true;
	omp_set_nested(1);
#endif

	if (argc != 4) {
		if (procID == 0) {
			std::cerr << "Bad arguments\nUsage: " << argv[0]
					<< " netFileName descFileName seed\n";
		}
		exit(1);
	}

	if (procID == 0) {
		std::cout << "#";
		for (int i = 0; i < argc; i++) {
			std::cout << argv[i] << " ";
		}
		std::cout << std::endl;
	}
	auto start = std::chrono::steady_clock::now();

	std::string netFileName(argv[1]);
	std::string fileName(argv[2]);
	long int seed = std::atol(argv[3]);
	bool keepConnected;
	double fnRatio;
	double tnRatio;
	std::size_t nbTests;
	double ratioStart;
	double ratioStep;
	double ratioEnd;

	std::vector<PerfM> perfMeasures;
	std::vector<std::string> predictors;

	readFile(fileName, keepConnected, fnRatio, tnRatio, nbTests, ratioStart,
			ratioStep, ratioEnd, perfMeasures, predictors, procID);

	RandomGen rng(seed);
	auto refNet = DNetwork<>::read(netFileName, false, true);
	std::size_t n = refNet->getNbNodes();
	std::size_t m = refNet->getNbEdges();
	if (procID == 0) {
		std::cout << "# n: " << n << " m: " << m << std::endl;
	}
	if (keepConnected) {
		ratioEnd = std::min(ratioEnd, (m - (n - 1.0)) / m);
	}
	bool first = true;
	for (std::size_t i = 0; i < nbTests; i++) {
		for (double ratio = ratioStart; ratio <= (ratioEnd + 1.0e-6); ratio +=
				ratioStep) {
			auto testData = NetworkManipulator<DNetwork<>>::createTestData(
					refNet, ratio, 0, keepConnected, false, fnRatio, 0,
					tnRatio, FN, TN, rng.getInt(), false);

			std::vector<
					std::shared_ptr<
							PerfMeasure<
									PredResults<TestData<DNetwork<>>,
											DLPredictor<>>>>> pms;
			//  performance measure
			for (auto it = perfMeasures.begin(); it != perfMeasures.end();
					++it) {
				auto pm = createPerfM(*it, testData.getNbPos()); //Assumes pos already generated
#ifdef WITH_OPENMP
						pm->setParallel(parallel);
#endif
				if (pm->requiresPos()) {
					testData.genPos();
				}
				if (pm->requiresNeg()) {
					testData.genNeg();
				}
				pms.push_back(pm);
			}

			testData.lock();

			PerfEvaluator<TestData<DNetwork<>>, DLPredictor<>> perf(testData);

			// predictors 
			for (auto it = predictors.begin(); it != predictors.end(); ++it) {
				auto predictor = createPredictor(*it, testData, rng);
#ifdef WITH_OPENMP
				predictor->setParallel(parallel);
#endif
#ifdef WITH_MPI
				predictor->setDistributed(distributed);
#endif
				perf.addPredictor(predictor);
			}

			//  performance measure
			for (auto it = pms.begin(); it != pms.end(); ++it) {
				perf.addPerfMeasure(*it);
			}
			perf.eval();

			if (procID == 0) {
				if (first) {
					std::cout << "#ratio\t";
					for (auto it = perf.resultsBegin(); it != perf.resultsEnd();
							++it) {
						std::cout << it->first << "\t";
					}
					std::cout << std::endl;
					first = false;
				}
				std::cout << std::fixed << std::setprecision(2) << ratio
						<< "\t";
				for (auto it = perf.resultsBegin(); it != perf.resultsEnd();
						++it) {
					std::cout << std::fixed << std::setprecision(4)
							<< it->second << "\t";
				}
				std::cout << std::endl;
			}
			refNet->shuffle(rng.getInt());
		}
	}

	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;
	if (procID == 0) {
		std::cerr << "#Time: "
				<< std::chrono::duration<double, std::milli>(diff).count()
				<< " ms" << std::endl;
	}

#ifdef WITH_MPI
	MPI_Finalize();
#endif

	return 0;
}

