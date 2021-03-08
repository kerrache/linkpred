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

#include "linkpred/simp/evaluator.hpp"

namespace LinkPred {
namespace Simp {

Evaluator::Evaluator() {
	factory = std::make_shared<Factory>();
}

std::shared_ptr<Encoder<>> Evaluator::Factory::createEncoder(std::string name,
		std::shared_ptr<const LinkPred::UNetwork<>> net, int dim,
		long int seed) {
	std::shared_ptr<Encoder<>> encoder;
	if (name == "DPW") {
		encoder = std::make_shared<DeepWalk<>>(net, seed);
		if (dim > 0) {
			encoder->setDim(dim);
		}
	} else if (name == "HMSM") {
		encoder = std::make_shared<HMSM<>>(net, seed);
		if (dim > 0) {
			encoder->setDim(dim);
		}
	} else if (name == "LVS") {
		encoder = std::make_shared<LargeVis<>>(net, seed);
		if (dim > 0) {
			encoder->setDim(dim);
		}
	} else
#ifdef LINKPRED_WITH_ARMADILLO
		if (name == "LEM") {
		encoder = std::make_shared<LEM<>>(net);
		if (dim > 0) {
			encoder->setDim(dim);
		}
	} else if (name == "LLE") {
		encoder = std::make_shared<LLE<>>(net);
		if (dim > 0) {
			encoder->setDim(dim);
		}
	} else
#endif
		if (name == "LIN") {
		encoder = std::make_shared<LINE<>>(net, seed);
		if (dim > 0) {
			encoder->setDim(dim);
		}
	} else if (name == "MFC") {
		encoder = std::make_shared<MatFact<>>(net, seed);
		if (dim > 0) {
			encoder->setDim(dim);
		}
	} else if (name == "N2V") {
		encoder = std::make_shared<Node2Vec<>>(net, seed);
		if (dim > 0) {
			encoder->setDim(dim);
		}
	} else {
		throw std::invalid_argument("Unknown encoder name");
	}
	return encoder;
}

std::shared_ptr<Classifier<>> Evaluator::Factory::createClassifier(std::string name,
		int dim, long int seed) {

	std::shared_ptr<Classifier<>> classifier;
#ifdef LINKPRED_WITH_MLPACK
	if (name == "FFN") {
		auto ffn = std::make_shared<FFN<>>();
		ffn->setAutoArch(dim);
		classifier = ffn;
	} else if (name == "LSVM") {
		classifier = std::make_shared<LinearSVM<>>();

	} else if (name == "NVB") {
		classifier = std::make_shared<NaiveBayes<>>();
	} else
#endif
	if (name == "LGR") {
		classifier = std::make_shared<LogisticRegresser<>>(0.001, seed);
	} else {
		throw std::invalid_argument("Unknown classifier name");
	}
	return classifier;
}

std::shared_ptr<SimMeasure> Evaluator::Factory::createSimMeasure(std::string name) {
	std::shared_ptr<SimMeasure> simMeasure;

	if (name == "CSM") {
		simMeasure = std::make_shared<CosineSim>();
	} else if (name == "DTP") {
		simMeasure = std::make_shared<DotProd>();
	} else if (name == "L1") {
		simMeasure = std::make_shared<L1Sim>();
	} else if (name == "L2") {
		simMeasure = std::make_shared<L2Sim>();
	} else if (name == "PRS") {
		simMeasure = std::make_shared<Pearson>();
	} else {
		throw std::invalid_argument("Unknown similarity measure name");
	}

	return simMeasure;
}

std::vector<std::shared_ptr<ULPredictor<>>> Evaluator::Factory::getPredictors(
		std::shared_ptr<UNetwork<> const> obsNet) {
	for (auto it = ada.begin(); it != ada.end(); ++it) {
		auto pr = std::make_shared<UADAPredictor<>>(obsNet);
		pr->setName(*it);
		predictors.push_back(pr);
	}

	for (auto it = cne.begin(); it != cne.end(); ++it) {
		auto pr = std::make_shared<UCNEPredictor<>>(obsNet);
		pr->setName(*it);
		predictors.push_back(pr);
	}

	for (auto it = cra.begin(); it != cra.end(); ++it) {
		auto pr = std::make_shared<UCRAPredictor<>>(obsNet);
		pr->setName(*it);
		predictors.push_back(pr);
	}

	for (auto it = ecl.begin(); it != ecl.end(); ++it) {

		auto encoderName = it->second.encoderName;
		auto classifierName = it->second.classifierName;
		auto dim = it->second.dim;
		auto posRatio = it->second.posRatio;
		auto negRatio = it->second.negRatio;
		auto seed = it->second.seed;

		RandomGen rng(seed);
		auto encoder = createEncoder(encoderName, obsNet, dim, rng.getInt());
		auto classifier = createClassifier(classifierName,
				encoder->getEdgeCodeDim(), rng.getInt());
		auto pr = std::make_shared<UECLPredictor<>>(obsNet, encoder,
				classifier, rng.getInt());

		pr->setPosRatio(posRatio);
		pr->setNegRatio(negRatio);
		pr->setName(it->first);

		predictors.push_back(pr);
	}

	for (auto it = esm.begin(); it != esm.end(); ++it) {

		auto encoderName = it->second.encoderName;
		auto simMeasureName = it->second.simMeasureName;
		auto dim = it->second.dim;
		auto seed = it->second.seed;

		RandomGen rng(seed);
		auto encoder = createEncoder(encoderName, obsNet, dim, rng.getInt());
		auto simMeasure = createSimMeasure(simMeasureName);
		auto pr = std::make_shared<UESMPredictor<>>(obsNet, encoder,
				simMeasure);
		pr->setName(it->first);

		predictors.push_back(pr);
	}

	for (auto it = fbm.begin(); it != fbm.end(); ++it) {
		auto pr = std::make_shared<UFBMPredictor<>>(obsNet,
				it->second.seed);
		pr->setName(it->first);
		pr->setMaxIter(it->second.maxIter);
		predictors.push_back(pr);
	}

	for (auto it = hdi.begin(); it != hdi.end(); ++it) {
		auto pr = std::make_shared<UHDIPredictor<>>(obsNet);
		pr->setName(*it);
		predictors.push_back(pr);
	}

	for (auto it = hpi.begin(); it != hpi.end(); ++it) {
		auto pr = std::make_shared<UHPIPredictor<>>(obsNet);
		pr->setName(*it);
		predictors.push_back(pr);
	}

	for (auto it = hrg.begin(); it != hrg.end(); ++it) {
		auto pr = std::make_shared<UHRGPredictor<>>(obsNet,
				it->second.seed);
		pr->setName(it->first);
		pr->setNbBeans(it->second.nbBeans);
		pr->setNbSamples(it->second.nbSamples);
		predictors.push_back(pr);
	}

	for (auto it = hyp.begin(); it != hyp.end(); ++it) {
		auto pr = std::make_shared<UHYPPredictor<>>(obsNet,
				it->second.seed);
		pr->setName(it->first);
		pr->setM(it->second.m);
		pr->setL(it->second.L);
		pr->setGamma(it->second.gamma);
		pr->setZeta(it->second.zeta);
		pr->setT(it->second.T);
		predictors.push_back(pr);
	}

	for (auto it = jid.begin(); it != jid.end(); ++it) {
		auto pr = std::make_shared<UJIDPredictor<>>(obsNet);
		pr->setName(*it);
		predictors.push_back(pr);
	}

	for (auto it = kab.begin(); it != kab.end(); ++it) {
		auto pr = std::make_shared<UKABPredictor<>>(obsNet);
		pr->setName(it->first);
		pr->setHorizLim(it->second.horizLim);
		predictors.push_back(pr);
	}

	for (auto it = lcp.begin(); it != lcp.end(); ++it) {
		auto pr = std::make_shared<ULCPPredictor<>>(obsNet);
		pr->setName(it->first);
		pr->setEpsilon(it->second.epsilon);
		predictors.push_back(pr);
	}

	for (auto it = lhn.begin(); it != lhn.end(); ++it) {
		auto pr = std::make_shared<ULHNPredictor<>>(obsNet);
		pr->setName(*it);
		predictors.push_back(pr);
	}

	for (auto it = pat.begin(); it != pat.end(); ++it) {
		auto pr = std::make_shared<UPATPredictor<>>(obsNet);
		pr->setName(*it);
		predictors.push_back(pr);
	}

	for (auto it = pst.begin(); it != pst.end(); ++it) {
		auto pr = std::make_shared<UPSTPredictor<>>(obsNet);
		pr->setName(it->first);
		pr->loadEdgeScores(it->second.fileName);
		predictors.push_back(pr);
	}

	for (auto it = ral.begin(); it != ral.end(); ++it) {
		auto pr = std::make_shared<URALPredictor<>>(obsNet);
		pr->setName(*it);
		predictors.push_back(pr);
	}

	for (auto it = rnd.begin(); it != rnd.end(); ++it) {
		auto pr = std::make_shared<URNDPredictor<>>(obsNet,
				it->second.seed);
		pr->setName(it->first);
		predictors.push_back(pr);
	}

	for (auto it = sai.begin(); it != sai.end(); ++it) {
		auto pr = std::make_shared<USAIPredictor<>>(obsNet);
		pr->setName(*it);
		predictors.push_back(pr);
	}

	for (auto it = sbm.begin(); it != sbm.end(); ++it) {
		auto pr = std::make_shared<USBMPredictor<>>(obsNet,
				it->second.seed);
		pr->setName(it->first);
		pr->setMaxIter(it->second.maxIter);
		predictors.push_back(pr);
	}

	for (auto it = shp.begin(); it != shp.end(); ++it) {
		auto pr = std::make_shared<USHPPredictor<>>(obsNet,
				it->second.seed);
		pr->setName(it->first);
		predictors.push_back(pr);
	}

	for (auto it = soi.begin(); it != soi.end(); ++it) {
		auto pr = std::make_shared<USOIPredictor<>>(obsNet);
		pr->setName(*it);
		predictors.push_back(pr);
	}

	if (predictors.size() == 0) {
		throw std::runtime_error("No predictors were added.");
	}

#ifdef LINKPRED_WITH_OPENMP
	for (std::size_t i = 0; i < predictors.size(); i++) {
		predictors[i]->setParallel(parallel);
	}
#endif

	return predictors;
}

std::vector<std::shared_ptr<PerfMeasure<>>> Evaluator::Factory::getPerfMeasures(
		TestData<> const &testData) {
	std::vector < std::shared_ptr<PerfMeasure<>> > pms;
	if (enableROC) {
		pms.push_back(std::make_shared<ROC<>>());
	}
	if (enablePR) {
		pms.push_back(std::make_shared<PR<>>());
	}
	if (enableTPR) {
		pms.push_back(std::make_shared<TPR<>>(testData.getNbPos()));
	}
	if (pms.size() == 0) {
		throw std::runtime_error("No performances measures were added.");
	}
	return pms;
}

void Evaluator::addADA(std::string const & name) {
	if (factory->ada.find(name) != factory->ada.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->ada.insert(name);
}

void Evaluator::addCNE(std::string const & name) {
	if (factory->cne.find(name) != factory->cne.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->cne.insert(name);
}

void Evaluator::addCRA(std::string const & name) {
	if (factory->cra.find(name) != factory->cra.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->cra.insert(name);
}

void Evaluator::addECL(std::string const & name, std::string encoderName,
		std::string classifierName, int dim, double posRatio, double negRatio,
		long int seed) {
	if (factory->ecl.find(name) != factory->ecl.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->ecl[name] = {encoderName, classifierName, dim, posRatio, negRatio,
		seed};
}

void Evaluator::addESM(std::string const & name, std::string encoderName,
		std::string simMeasureName, int dim, long int seed) {
	if (factory->esm.find(name) != factory->esm.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->esm[name] = {encoderName, simMeasureName, dim, seed};
}

void Evaluator::addFBM(std::string const & name, int maxIter, long int seed) {
	if (factory->fbm.find(name) != factory->fbm.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->fbm[name] = {maxIter, seed};
}

void Evaluator::addHDI(std::string const & name) {
	if (factory->hdi.find(name) != factory->hdi.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->hdi.insert(name);
}

void Evaluator::addHPI(std::string const & name) {
	if (factory->hpi.find(name) != factory->hpi.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->hpi.insert(name);
}

void Evaluator::addHRG(std::string const & name, int nbBeans, int nbSamples,
		long int seed) {
	if (factory->hrg.find(name) != factory->hrg.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->hrg[name] = {nbBeans, nbSamples, seed};
}

void Evaluator::addHYP(std::string const & name, double m, double L,
		double gamma, double zeta, double T, long int seed) {
	if (factory->hyp.find(name) != factory->hyp.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->hyp[name] = {m, L, gamma, zeta, T, seed};
}

void Evaluator::addJID(std::string const & name) {
	if (factory->jid.find(name) != factory->jid.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->jid.insert(name);
}

void Evaluator::addKAB(std::string const & name, int horizLim) {
	if (factory->kab.find(name) != factory->kab.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->kab[name] = {horizLim};
}

void Evaluator::addLCP(std::string const & name, double epsilon) {
	if (factory->lcp.find(name) != factory->lcp.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->lcp[name] = {epsilon};
}

void Evaluator::addLHN(std::string const & name) {
	if (factory->lhn.find(name) != factory->lhn.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->lhn.insert(name);
}

void Evaluator::addPAT(std::string const & name) {
	if (factory->pat.find(name) != factory->pat.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->pat.insert(name);
}

void Evaluator::addPST(std::string const & name, std::string fileName) {
	if (factory->pst.find(name) != factory->pst.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->pst[name] = {fileName};
}

void Evaluator::addRAL(std::string const & name) {
	if (factory->ral.find(name) != factory->ral.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->ral.insert(name);
}

void Evaluator::addRND(std::string const & name, long int seed) {
	if (factory->rnd.find(name) != factory->rnd.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->rnd[name] = {seed};
}

void Evaluator::addSAI(std::string const & name) {
	if (factory->sai.find(name) != factory->sai.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->sai.insert(name);
}

void Evaluator::addSBM(std::string const & name, int maxIter, long int seed) {
	if (factory->sbm.find(name) != factory->sbm.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->sbm[name] = {maxIter, seed};
}

void Evaluator::addSHP(std::string const & name, long int seed) {
	if (factory->shp.find(name) != factory->shp.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->shp[name] = {seed};
}

void Evaluator::addSOI(std::string const & name) {
	if (factory->soi.find(name) != factory->soi.end()) {
		throw std::runtime_error("Predictor " + name + " already exists");
	}
	factory->soi.insert(name);
}

void Evaluator::addROC() {
	factory->enableROC = true;
}

void Evaluator::addPR() {
	factory->enablePR = true;
}

void Evaluator::addTPR() {
	factory->enableTPR = true;
}

void Evaluator::genTestData(std::string const & fullNetFileName,
		std::string const &obsEdgesFileName,
		std::string const &remEdgesFileName, double remRatio,
		bool keepConnected, long int seed) {

	auto refNet = UNetwork<>::read(fullNetFileName);

	auto testData = NetworkManipulator<>::createTestData(refNet, remRatio,
			0, keepConnected, false, 1, 0, 1, LinkPred::FN,
			LinkPred::TN, seed, true);

	{
		std::ofstream out;
		out.open(obsEdgesFileName.c_str(), std::fstream::out);

		if (!out) {
			throw std::runtime_error("Cannot open file: " + obsEdgesFileName);
		}

		auto obsNet = testData.getObsNet();
		for (auto it = obsNet->edgesBegin(); it != obsNet->edgesEnd(); ++it) {
			auto i = obsNet->getLabel(obsNet->start(*it));
			auto j = obsNet->getLabel(obsNet->end(*it));
			out << i << "\t" << j << std::endl;
		}

		out.close();
	}

	{
		std::ofstream out;
		out.open(remEdgesFileName.c_str(), std::fstream::out);

		if (!out) {
			throw std::runtime_error("Cannot open file: " + remEdgesFileName);
		}

		auto obsNet = testData.getObsNet();
		for (auto it = testData.posBegin(); it != testData.posEnd(); ++it) {
			auto i = obsNet->getLabel(obsNet->start(*it));
			auto j = obsNet->getLabel(obsNet->end(*it));
			out << i << "\t" << j << std::endl;
		}

		out.close();
	}
}

void Evaluator::run(std::string const &fullNetFileName, int nbRuns,
		double remRatio, bool keepConnected, long int seed) {

	auto net = UNetwork<>::read(fullNetFileName, false, true);
	PerfeEvalExpDescp<> ped;
	ped.refNet = net;
	ped.keepConnected = keepConnected;
	ped.nbTestRuns = nbRuns;
	ped.seed = seed;
	PerfEvalExp<> exp(ped, factory);
	exp.run();

	perfRes.clear();
	perfRes.reserve(nbRuns);
	for (auto it = exp.resultsBegin(); it != exp.resultsEnd(); ++it) {
		std::vector<PerfRes> vals;
		for (auto rit = it->begin(); rit != it->end(); ++rit) {
			vals.push_back( { rit->first, rit->second });
		}
		perfRes.push_back(std::move(vals));
	}
}

void Evaluator::run(std::string const &obsEdgesFileName,
		std::string const &remEdgesFileName) {

	auto testData = NetworkManipulator<>::loadTestData(obsEdgesFileName,
			remEdgesFileName, "", true, 1, true, 1, LinkPred::FN,
			LinkPred::TN, 0, false);
	auto predictors = factory->getPredictors(testData.getObsNet());
	auto perfMeasures = factory->getPerfMeasures(testData);
	for (std::size_t i = 0; i < perfMeasures.size(); i++) {
		if (perfMeasures[i]->requiresPos()) {
			testData.genPos();
		}
		if (perfMeasures[i]->requiresNeg()) {
			testData.genNeg();
		}
	}
	testData.lock();
	PerfEvaluator<> perf(testData);

	for (std::size_t i = 0; i < predictors.size(); i++) {
		perf.addPredictor(predictors[i]);
	}

	for (std::size_t i = 0; i < perfMeasures.size(); i++) {
		perf.addPerfMeasure(perfMeasures[i]);
	}

	perf.eval();

	perfRes.clear();
	perfRes.reserve(1);
	std::vector<PerfRes> vals;
	for (auto it = perf.resultsBegin(); it != perf.resultsEnd(); ++it) {
		vals.push_back( { it->first, it->second });
	}

	for (auto it = vals.begin(); it != vals.end(); ++it) {
		std::cout << it->name << "\t";
	}
	std::cout << std::endl;
	for (auto it = vals.begin(); it != vals.end(); ++it) {
		std::cout << std::fixed << std::setprecision(4) << it->res << "\t";
	}
	std::cout << std::endl;

	perfRes.push_back(std::move(vals));
}

}
/* namespace Simp */
}
/* namespace LinkPred */

