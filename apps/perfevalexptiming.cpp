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
	ADA, /**< Adamic Adar predictor. */
	CNE, /**< Common neighbors. */
	CRA, /**< CRA. */
	CST, /**< Constant. */
	FBM, /**< Fast blocking model. */
	HDI, /**< Hub depromoted index predictor. */
	HPI, /**< Hub promoted index predictor. */
	HRG, /**< HRG. */
	HYP, /**< Hypermap. */
	JID, /**< Jackard index. */
	LCP, /**< Local path predictor. */
	LHN, /**< Leicht-Holme-Newman index predictor. */
	PAT, /**< Preferential attachment predictor. */
	RAL, /**< Resource allocation predictor. */
	RND, /**< Random. */
	SAI, /**< Salton index predictor. */
	SBM, /**< Stochastic bloc model. */
	SHP, /**< Shortest path. */
	SOI, /**< Sorensen index predictor. */
	KAB, /**< Scalable popularity-similarity. */
	SUM, /**< Sum of degrees. */
	POP, /**< Product and sum of degrees. */
	NED, /**< Degrees of neighbors. */
	PND  /**< Popularity-degrees of neighbors. */
};

std::map<std::string, LPMethod> lpNameMap = { { "ADA", ADA }, { "CNE", CNE }, { "CRA", CRA }, { "CST", CST }, { "FBM", FBM }, { "HDI", HDI }, { "HPI", HPI }, { "HRG", HRG }, { "HYP", HYP }, { "JID", JID }, { "LCP", LCP }, { "LHN", LHN }, { "PAT", PAT }, { "RAL", RAL }, { "RND", RND }, { "SAI", SAI }, { "SBM", SBM }, { "SHP", SHP }, { "SOI", SOI }, { "KAB", KAB }, { "SUM", SUM }, { "POP", POP }, { "NED", NED }, { "PND", PND } };

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

std::shared_ptr<ULPredictor<>> createPredictor(std::string lpName,
		std::shared_ptr<UNetwork<> const> obsNet, RandomGen & rng) {

	auto tokens = split(lpName, '_');
	auto lpm = lpNameMap.at(tokens[0]);
	switch (lpm) {

	case CNE:
		return std::make_shared<UCNEPredictor<>>(obsNet); /**< Common neighbors. */

	case CRA:
		return std::make_shared<UCRAPredictor<>>(obsNet); /**< CRA. */

	case CST:
		return std::make_shared<UCSTPredictor<>>(obsNet, rng.getInt()); /**< CST. */

	case JID:
		return std::make_shared<UJIDPredictor<>>(obsNet); /**< Jackard index. */

	case SBM:
		return std::make_shared<USBMPredictor<>>(obsNet, rng.getInt()); /**< Stochastic bloc model. */

	case SHP:
		return std::make_shared<USHPPredictor<>>(obsNet, rng.getInt()); /**< Stochastic bloc model. */

	case HRG:
		return std::make_shared<UHRGPredictor<>>(obsNet, rng.getInt()); /**< HRG. */

	case KAB:
		return std::make_shared<UKABPredictor<>>(obsNet, rng.getInt()); /**< Scalable popularity-similarity. */

	case FBM:
		return std::make_shared<UFBMPredictor<>>(obsNet, rng.getInt()); /**< Fast blocking model. */

	case HYP:
		return std::make_shared<UHYPPredictor<>>(obsNet, rng.getInt()); /**< Hypermap. */

	case ADA:
		return std::make_shared<UADAPredictor<>>(obsNet); /**< Adamic Adar solver. */

	case RAL:
		return std::make_shared<URALPredictor<>>(obsNet); /**< Resource allocation solver. */

	case PAT:
		return std::make_shared<UPATPredictor<>>(obsNet); /**< Preferential attachment solver. */

	case LCP:
		return std::make_shared<ULCPPredictor<>>(obsNet); /**< Local path solver. */

	case RND:
		return std::make_shared<URNDPredictor<>>(obsNet, rng.getInt()); /**< Random solver. */

	case SAI:
		return std::make_shared<USAIPredictor<>>(obsNet); /**< Salton index solver. */

	case SOI:
		return std::make_shared<USOIPredictor<>>(obsNet); /**< Sorensen index solver. */

	case HPI:
		return std::make_shared<UHPIPredictor<>>(obsNet); /**< Hub promoted index solver. */

	case HDI:
		return std::make_shared<UHDIPredictor<>>(obsNet); /**< Hub depromoted index solver. */

	case LHN:
		return std::make_shared<ULHNPredictor<>>(obsNet); /**< Leicht-Holme-Newman index solver. */

	case SUM:
		return std::make_shared<USUMPredictor<>>(obsNet);

	case POP:
		return std::make_shared<UPOPPredictor<>>(obsNet);

	case NED:
		return std::make_shared<UNEDPredictor<>>(obsNet);

	case PND:
		return std::make_shared<UPNDPredictor<>>(obsNet);

	default:
		throw std::invalid_argument("Unknown link prediction method");
	}
}

void readDescFile(std::string fileName, PerfeEvalExpDescp<> & ped,
		std::vector<PerfM> & perfMeasures,
		std::vector<std::string> & predictors, int procID) {

	bool keepConnected;
	double fnRatio;
	double tnRatio;
	std::size_t nbTests;
	double ratioStart;
	double ratioStep;
	double ratioEnd;

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

	ped.nbTestRuns = nbTests; /**< Number of test runs. */
	ped.keepConnected = keepConnected; /**< Whether to keep the network connected. */
	ped.fnRatio = fnRatio; /**< Ratio of false negatives used in the test set. */
	ped.tnRatio = tnRatio; /**< Ratio of true negatives used in the test set. */
	ped.ratioStart = ratioStart; /**< Start value of the ratio of removed edges. */
	ped.ratioEnd = ratioEnd; /**< End value of the ratio of removed edges. This is adjusted if keepConnected is set to true. */
	ped.ratioStep = ratioStep; /**< Step size of the ratio of removed edges. */
}

std::shared_ptr<PerfMeasure<>> createPerfM(PerfM pm, long int nb) {
	switch (pm) {
	case ROCT:
		return std::make_shared<ROC<>>();
	case PRT:
		return std::make_shared<PR<>>();
	case TPRT:
		return std::make_shared<TPR<>>(nb, "TPR");
	case TPRTT: {
		auto tpr = std::make_shared<TPR<>>(nb, "TPRT");
		tpr->setUseTopMethod(true);
		return tpr;
	}
	default:
		throw std::invalid_argument("Unknown performance measure");
	}
}

class Factory: public PEFactory<> {
protected:
	std::vector<PerfM> perfMeasures;
	std::vector<std::string> predictors;
	bool parallel = false;
	bool distributed = true;
	RandomGen rng;
public:
	Factory(std::vector<PerfM> const & perfMeasures,
			std::vector<std::string> const &predictors, bool parallel,
			bool distributed, long int seed) :
			perfMeasures(perfMeasures), predictors(predictors), parallel(
					parallel), distributed(distributed), rng(seed) {
	}

	virtual std::vector<std::shared_ptr<ULPredictor<>>> getPredictors(
			std::shared_ptr<UNetwork<> const> obsNet) {
		std::vector<std::shared_ptr<ULPredictor<>>> prs;
		for (auto it = predictors.begin(); it != predictors.end(); ++it) {
			auto predictor = createPredictor(*it, obsNet, rng);
#ifdef WITH_OPENMP
			predictor->setParallel(parallel);
#endif
#ifdef WITH_MPI
			predictor->setDistributed(distributed);
#endif
			prs.push_back(predictor);
		}
		return prs;
	}

	virtual std::vector<std::shared_ptr<PerfMeasure<>>> getPerfMeasures(
			TestData<> const & testData) {
		std::vector<std::shared_ptr<PerfMeasure<>>> pms;
		for (auto it = perfMeasures.begin(); it != perfMeasures.end(); ++it) {
			auto pm = createPerfM(*it, testData.getNbPos()); //Assumes pos already generated
#ifdef WITH_OPENMP
			pm->setParallel(parallel);
#endif
			pms.push_back(pm);
		}
		return pms;
	}

	/**
	 * Destructor.
	 */
	virtual ~Factory() = default;
};

int main(int argc, char*argv[]) {

	int procID = 0;
#ifdef WITH_MPI
	bool distributed = true;
	int nbProcs = 1;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nbProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);
#else
	bool distributed = false;
#endif
#ifdef WITH_OPENMP
	bool parallel = true;
	omp_set_nested(1);
#else
	bool parallel = false;
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

	std::string netFileName(argv[1]);
	auto refNet = UNetwork<>::read(netFileName, false, true);
	std::string descFileName(argv[2]);
	PerfeEvalExpDescp<> ped;
	ped.refNet = refNet; /**< Reference network. */
	std::vector<PerfM> perfMeasures;
	std::vector<std::string> predictors;
	readDescFile(descFileName, ped, perfMeasures, predictors, procID);
	ped.timingEnabled = true; // Enable timing
	long int seed = std::atol(argv[3]);
	RandomGen rng(seed);
	ped.seed = rng.getInt(); /**< Seed for the random number generator. */
#ifdef WITH_MPI
	ped.distributed = distributed;
#endif
	auto factory = std::make_shared<Factory>(perfMeasures, predictors, parallel,
			distributed, rng.getInt());
	PerfEvalExp<> exp(ped, factory);
	exp.run();
#ifdef WITH_MPI
	MPI_Finalize();
#endif

	return 0;
}

