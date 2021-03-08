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

#include <linkpred/utils/miscutils.hpp>
#include "linkpred/perf/perfevaluator.hpp"
#include <iostream>
#include <limits>

namespace LinkPred {


template<typename TestDataT, typename LPredictorT, typename PredResultsT,
		typename PerfMeasureT> void PerfEvaluator<TestDataT, LPredictorT,
		PredResultsT, PerfMeasureT>::evalNoTiming() {
	logger(logDebug, "Computing performance measures without timing...")

	std::vector<PerfResults> tmpRes;
	tmpRes.resize(predictors.size());
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (std::size_t i = 0; i < predictors.size(); i++) {
		auto predictor = predictors[i];
		predictor->init();
		predictor->learn();
		auto predResults = std::make_shared<PredResultsT>(testData, predictor);
		for (std::size_t j = 0; j < measures.size(); j++) {
			auto measure = measures[j];
			PerfResults res;
			measure->eval(predResults, res);
			for (auto it = res.begin(); it != res.end(); ++it) {
				tmpRes[i][it->first + predictor->getName()] = it->second;
			}
		}

		// Free some memory if necessary
		predictors[i].reset();
	}
	for (std::size_t i = 0; i < predictors.size(); i++) {
		for (auto it = tmpRes[i].begin(); it != tmpRes[i].end(); ++it) {
			results[it->first] = it->second;
		}
	}
	logger(logDebug, "Done")

}

template<typename TestDataT, typename LPredictorT, typename PredResultsT,
		typename PerfMeasureT> void PerfEvaluator<TestDataT, LPredictorT,
		PredResultsT, PerfMeasureT>::evalTiming() {
	logger(logDebug, "Computing performance measures with timing...")

	std::vector<PerfResults> tmpRes;
	tmpRes.resize(predictors.size());
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (std::size_t i = 0; i < predictors.size(); i++) {

		auto predictor = predictors[i];

		// Initializing predictors
		{
			auto start = std::chrono::steady_clock::now();
			predictor->init();
			auto end = std::chrono::steady_clock::now();

			auto diff = end - start;
			tmpRes[i]["ITN" + predictor->getName()] = std::chrono::duration<
					double, std::nano>(diff).count();
		}

		// Learning
		{
			auto start = std::chrono::steady_clock::now();
			predictor->learn();
			auto end = std::chrono::steady_clock::now();

			auto diff = end - start;
			tmpRes[i]["LTN" + predictor->getName()] = std::chrono::duration<
					double, std::nano>(diff).count();
		}

		auto predResults = std::make_shared<PredResultsT>(testData, predictor);

		// Predicting
		{
			auto start = std::chrono::steady_clock::now();
			predResults->compPosScores();
			predResults->compNegScores();
			auto end = std::chrono::steady_clock::now();

			auto diff = end - start;
			tmpRes[i]["PTN" + predictor->getName()] = std::chrono::duration<
					double, std::nano>(diff).count();
		}

		tmpRes[i]["TTN" + predictor->getName()] = tmpRes[i]["ITN"
				+ predictor->getName()]
				+ tmpRes[i]["LTN" + predictor->getName()]
				+ tmpRes[i]["PTN" + predictor->getName()];

		for (std::size_t j = 0; j < measures.size(); j++) {
			auto measure = measures[j];
			PerfResults res;
			measure->eval(predResults, res);
			for (auto it = res.begin(); it != res.end(); ++it) {
				tmpRes[i][it->first + predictor->getName()] = it->second;
			}
		}

		// Free some memory if necessary
		predictors[i].reset();
	}
	for (std::size_t i = 0; i < predictors.size(); i++) {
		for (auto it = tmpRes[i].begin(); it != tmpRes[i].end(); ++it) {
			results[it->first] = it->second;
		}
	}
	logger(logDebug, "Done")
}

template<typename Network, typename TestDataT, typename LPredictorT,
		typename PredResultsT, typename PerfMeasureT> void PerfEvalExp<Network,
		TestDataT, LPredictorT, PredResultsT, PerfMeasureT>::run() {
	logger(logDebug, "Computing performance measures...")
	if (ped.timingEnabled) {
		runTiming();
	} else {
		runNoTiming();
	}
	logger(logDebug, "Done")
}

template<typename Network, typename TestDataT, typename LPredictorT,
		typename PredResultsT, typename PerfMeasureT> void PerfEvalExp<Network,
		TestDataT, LPredictorT, PredResultsT, PerfMeasureT>::runNoTiming() {
	logger(logDebug, "Computing performance measures without timing...")

	auto refNet = ped.refNet;
	bool keepConnected = ped.keepConnected;
	double fnRatio = ped.fnRatio;
	double tnRatio = ped.tnRatio;
	double ratioStart = ped.ratioStart;
	double ratioEnd = ped.ratioEnd;
	double ratioStep = ped.ratioStep;
	auto nbTestRuns = ped.nbTestRuns;
	auto out = ped.out;
#ifdef LINKPRED_WITH_OPENMP
	bool parTestRuns = ped.parTestRuns;
	bool parPredictors = ped.parPredictors;
#endif
#ifdef LINKPRED_WITH_MPI
	bool distributed = ped.distributed;
	auto comm = ped.comm;
#endif
	RandomGen rng(ped.seed);
	std::size_t n = refNet->getNbNodes();
	std::size_t m = refNet->getNbEdges();
//The size of the range must be equal to the number of threads
	int procID = 0;
#ifdef LINKPRED_WITH_MPI
	int nbProcs = 1;
	if (distributed) {
		MPI_Comm_size(comm, &nbProcs);
		MPI_Comm_rank(comm, &procID);
	}
#endif
	if (procID == 0) {
		std::cout << "# n: " << n << " m: " << m << std::endl;
	}
	if (keepConnected) {
		ratioEnd = std::min(ratioEnd, (m - (n - 1.0)) / m);
	}
	results.clear();
	results.resize(nbTestRuns);
	bool first = true;
	auto start = std::chrono::steady_clock::now();
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parTestRuns)
#endif
	for (std::size_t k = 0; k < nbTestRuns; k++) {
		for (double ratio = ratioStart;
				ratio <= (ratioEnd + std::numeric_limits<double>::epsilon());
				ratio += ratioStep) {
			auto testData = NetworkManipulator<>::createTestData(refNet, ratio,
					0, keepConnected, false, fnRatio, 0, tnRatio, FN, TN,
					rng.getInt(), false);

			//  performance measure
			auto measures = factory->getPerfMeasures(testData);
			bool shuffle = false;
			for (auto it = measures.begin(); it != measures.end(); ++it) {
				auto pm = *it;
				if (pm->requiresPos()) {
					testData.genPos();
				}
				if (pm->requiresNeg()) {
					testData.genNeg();
				}
				if (pm->requiresShuffling()) {
					shuffle = true;
				}
			}

			testData.lock();

			// predictors
			auto predictors = factory->getPredictors(testData.getObsNet());

			std::vector<PerfResults> tmpRes;
			tmpRes.resize(predictors.size());
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parPredictors)
#endif
			for (std::size_t i = 0; i < predictors.size(); i++) {
				auto predictor = predictors[i];
				predictor->init();
				predictor->learn();
				auto predResults = std::make_shared<PredResultsT>(testData,
						predictor);
				for (std::size_t j = 0; j < measures.size(); j++) {
					auto measure = measures[j];
					PerfResults res;
					measure->eval(predResults, res);
					for (auto it = res.begin(); it != res.end(); ++it) {
						tmpRes[i][it->first + predictor->getName()] =
								it->second;
					}
				}
				// Free some memory if necessary
				predictors[i].reset();
			}
			for (std::size_t i = 0; i < predictors.size(); i++) {
				for (auto it = tmpRes[i].begin(); it != tmpRes[i].end(); ++it) {
					results[k][it->first] = it->second;
				}
			}
			if (procID == 0) {
				if (first) {
					(*out) << "#ratio\t";
					for (auto it = results[k].begin(); it != results[k].end();
							++it) {
						(*out) << it->first << "\t";
					}
					(*out) << std::endl;
					first = false;
				}
				(*out) << std::fixed << std::setprecision(2) << ratio << "\t";
				for (auto it = results[k].begin(); it != results[k].end();
						++it) {
					(*out) << std::fixed << std::setprecision(outPrec)
							<< it->second << "\t";
				}
				(*out) << std::endl;
			}
			if (shuffle) {
				refNet->shuffle(rng.getInt());
			}
		}
	}

	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;
	if (procID == 0) {
		(*out) << "#Time: "
				<< std::chrono::duration<double, std::milli>(diff).count()
				<< " ms" << std::endl;
	}

	logger(logDebug, "Done")
}

template<typename Network, typename TestDataT, typename LPredictorT,
		typename PredResultsT, typename PerfMeasureT> void PerfEvalExp<Network,
		TestDataT, LPredictorT, PredResultsT, PerfMeasureT>::runTiming() {
	logger(logDebug, "Computing performance measures with timing...")

	auto refNet = ped.refNet;
	bool keepConnected = ped.keepConnected;
	double fnRatio = ped.fnRatio;
	double tnRatio = ped.tnRatio;
	double ratioStart = ped.ratioStart;
	double ratioEnd = ped.ratioEnd;
	double ratioStep = ped.ratioStep;
	auto nbTestRuns = ped.nbTestRuns;
	auto out = ped.out;
#ifdef LINKPRED_WITH_OPENMP
	bool parTestRuns = ped.parTestRuns;
	bool parPredictors = ped.parPredictors;
#endif
#ifdef LINKPRED_WITH_MPI
	bool distributed = ped.distributed;
	auto comm = ped.comm;
#endif
	RandomGen rng(ped.seed);
	std::size_t n = refNet->getNbNodes();
	std::size_t m = refNet->getNbEdges();
	//The size of the range must be equal to the number of threads
	int procID = 0;
#ifdef LINKPRED_WITH_MPI
	int nbProcs = 1;
	if (distributed) {
		MPI_Comm_size(comm, &nbProcs);
		MPI_Comm_rank(comm, &procID);
	}
#endif
	if (procID == 0) {
		std::cout << "# n: " << n << " m: " << m << std::endl;
	}
	if (keepConnected) {
		ratioEnd = std::min(ratioEnd, (m - (n - 1.0)) / m);
	}
	results.clear();
	results.resize(nbTestRuns);
	bool first = true;
	auto start = std::chrono::steady_clock::now();
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parTestRuns)
#endif
	for (std::size_t k = 0; k < nbTestRuns; k++) {
		for (double ratio = ratioStart;
				ratio <= (ratioEnd + std::numeric_limits<double>::epsilon());
				ratio += ratioStep) {
			auto testData = NetworkManipulator<>::createTestData(refNet, ratio,
					0, keepConnected, false, fnRatio, 0, tnRatio, FN, TN,
					rng.getInt(), false);

			//  performance measure
			bool shuffle = false;
			auto measures = factory->getPerfMeasures(testData);
			for (auto it = measures.begin(); it != measures.end(); ++it) {
				auto pm = *it;
				if (pm->requiresPos()) {
					testData.genPos();
				}
				if (pm->requiresNeg()) {
					testData.genNeg();
				}
				if (pm->requiresShuffling()) {
					shuffle = true;
				}
			}

			testData.lock();

			// predictors
			auto predictors = factory->getPredictors(testData.getObsNet());
			std::vector<PerfResults> tmpRes;
			tmpRes.resize(predictors.size());
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parPredictors)
#endif
			for (std::size_t i = 0; i < predictors.size(); i++) {
				auto predictor = predictors[i];
				// Initializing predictors
				{
					auto start = std::chrono::steady_clock::now();
					predictor->init();
					auto end = std::chrono::steady_clock::now();

					auto diff = end - start;
					tmpRes[i]["ITN" + predictor->getName()] =
							std::chrono::duration<double, std::nano>(diff).count();
				}

				// Learning
				{
					auto start = std::chrono::steady_clock::now();
					predictor->learn();
					auto end = std::chrono::steady_clock::now();

					auto diff = end - start;
					tmpRes[i]["LTN" + predictor->getName()] =
							std::chrono::duration<double, std::nano>(diff).count();
				}

				auto predResults = std::make_shared<PredResultsT>(testData,
						predictor);

				// Predicting
				{
					auto start = std::chrono::steady_clock::now();
					predResults->compPosScores();
					predResults->compNegScores();
					auto end = std::chrono::steady_clock::now();

					auto diff = end - start;
					tmpRes[i]["PTN" + predictor->getName()] =
							std::chrono::duration<double, std::nano>(diff).count();
				}

				tmpRes[i]["TTN" + predictor->getName()] = tmpRes[i]["ITN"
						+ predictor->getName()]
						+ tmpRes[i]["LTN" + predictor->getName()]
						+ tmpRes[i]["PTN" + predictor->getName()];

				for (std::size_t j = 0; j < measures.size(); j++) {
					auto measure = measures[j];
					PerfResults res;
					measure->eval(predResults, res);
					for (auto it = res.begin(); it != res.end(); ++it) {
						tmpRes[i][it->first + predictor->getName()] =
								it->second;
					}
				}

				// Free some memory if necessary
				predictors[i].reset();
			}

			for (std::size_t i = 0; i < predictors.size(); i++) {
				for (auto it = tmpRes[i].begin(); it != tmpRes[i].end(); ++it) {
					results[k][it->first] = it->second;
				}
			}

			if (procID == 0) {
				if (first) {
					(*out) << "#ratio\t";
					for (auto it = results[k].begin(); it != results[k].end();
							++it) {
						(*out) << it->first << "\t";
					}
					(*out) << std::endl;
					first = false;
				}
				(*out) << std::fixed << std::setprecision(2) << ratio << "\t";
				for (auto it = results[k].begin(); it != results[k].end();
						++it) {
					(*out) << std::fixed << std::setprecision(4) << it->second
							<< "\t";
				}
				(*out) << std::endl;
			}
			if (shuffle) {
				refNet->shuffle(rng.getInt());
			}
		}
	}

	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;
	if (procID == 0) {
		(*out) << "#Time: "
				<< std::chrono::duration<double, std::milli>(diff).count()
				<< " ms" << std::endl;
	}

	logger(logDebug, "Done")

}

#define PERFEVALUATOR_CPP
#include "linkpred/instantiations.hpp"
#undef PERFEVALUATOR_CPP


} /* namespace LinkPred */
