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

#include "linkpred/simp/predictor.hpp"

namespace LinkPred {
namespace Simp {

std::vector<EdgeScore> Predictor::predAll(
		std::shared_ptr<ULPredictor<>> predictor) {

#ifdef LINKPRED_WITH_OPENMP
	predictor->setParallel(parallel);
#endif
	predictor->init();
	predictor->learn();

	std::vector<double> scores(net->getNbNonEdges());
	auto range = predictor->predictNeg(scores.begin());
	std::vector<EdgeScore> edgeScores(net->getNbNonEdges());
	std::size_t k = 0;
	for (auto it = range.first; it != range.second; ++it, k++) {
		EdgeScore es;
		es.i = net->getLabel(net->start(*it));
		es.j = net->getLabel(net->end(*it));
		es.score = scores[k];
		edgeScores[k] = es;
	}
	return edgeScores;
}

std::vector<EdgeScoreByID> Predictor::predAllByID(
		std::shared_ptr<ULPredictor<>> predictor) {

#ifdef LINKPRED_WITH_OPENMP
	predictor->setParallel(parallel);
#endif
	predictor->init();
	predictor->learn();

	std::vector<double> scores(net->getNbNonEdges());
	auto range = predictor->predictNeg(scores.begin());
	std::vector<EdgeScoreByID> edgeScores(net->getNbNonEdges());
	std::size_t k = 0;
	for (auto it = range.first; it != range.second; ++it, k++) {
		EdgeScoreByID es;
		es.i = net->start(*it);
		es.j = net->end(*it);
		es.score = scores[k];
		edgeScores[k] = es;
	}

	return edgeScores;
}

void Predictor::pred(std::shared_ptr<ULPredictor<>> predictor,
		std::vector<EdgeScore> &edgeScores) {

#ifdef LINKPRED_WITH_OPENMP
	predictor->setParallel(parallel);
#endif
	predictor->init();
	predictor->learn();

	std::vector<typename UNetwork<>::Edge> edges(edgeScores.size());
	std::vector<double> scores(edgeScores.size());
	std::size_t k = 0;
	for (auto it = edgeScores.begin(); it != edgeScores.end(); ++it, k++) {
		edges[k] = net->makeEdge(net->getID(it->i), net->getID(it->j));
	}
	predictor->predict(edges.begin(), edges.end(), scores.begin());
	for (k = 0; k < scores.size(); k++) {
		edgeScores[k].score = scores[k];
	}
}

std::vector<EdgeScore> Predictor::predTop(
		std::shared_ptr<ULPredictor<>> predictor, int k) {

#ifdef LINKPRED_WITH_OPENMP
	predictor->setParallel(parallel);
#endif
	predictor->init();
	predictor->learn();

	std::vector<typename UNetwork<>::Edge> edges(k);
	std::vector<double> scores(k);
	std::size_t nb = predictor->top(k, edges.begin(), scores.begin());
	std::vector<EdgeScore> edgeScores(nb);
	for (std::size_t i = 0; i < nb; i++) {
		EdgeScore es;
		es.i = net->getLabel(net->start(edges[i]));
		es.j = net->getLabel(net->end(edges[i]));
		es.score = scores[i];
		edgeScores[i] = es;
	}
	return edgeScores;
}

void Predictor::loadnet(std::string fileName) {
	net = UNetwork<>::read(fileName, false, true);
}

int Predictor::getNbNodes() const {
	return net->getNbNodes();
}

std::string Predictor::getLabel(int i) const {
	return net->getLabel(i);
}

int Predictor::getID(std::string const &i) const {
	return net->getID(i);
}

bool Predictor::isEdgeByID(int i, int j) const {
	return net->isEdge(i, j);
}

bool Predictor::isEdgeByLabel(std::string const &i, std::string const &j) const {
	return net->isEdge(net->getID(i), net->getID(j));
}

std::shared_ptr<Encoder<>> Predictor::createEncoder(std::string name,
		std::shared_ptr<const UNetwork<>> net, int dim, long int seed) {
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

std::shared_ptr<Classifier<>> Predictor::createClassifier(std::string name,
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

std::shared_ptr<SimMeasure> Predictor::createSimMeasure(std::string name) {
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

std::vector<EdgeScore> Predictor::predAllADA() {
	auto pr = std::make_shared<UADAPredictor<>>(net);
	return predAll(pr);
}

void Predictor::predADA(std::vector<EdgeScore> &edgeScores) {
	auto pr = std::make_shared<UADAPredictor<>>(net);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopADA(int k) {
	auto pr = std::make_shared<UADAPredictor<>>(net);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllCNE() {
	auto pr = std::make_shared<UCNEPredictor<>>(net);
	return predAll(pr);
}

void Predictor::predCNE(std::vector<EdgeScore> &edgeScores) {
	auto pr = std::make_shared<UCNEPredictor<>>(net);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopCNE(int k) {
	auto pr = std::make_shared<UCNEPredictor<>>(net);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllCRA() {
	auto pr = std::make_shared<UCRAPredictor<>>(net);
	return predAll(pr);
}

void Predictor::predCRA(std::vector<EdgeScore> &edgeScores) {
	auto pr = std::make_shared<UCRAPredictor<>>(net);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopCRA(int k) {
	auto pr = std::make_shared<UCRAPredictor<>>(net);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllECL(std::string encoderName,
		std::string classifierName, int dim, double posRatio, double negRatio,
		long int seed) {

	RandomGen rng(seed);
	auto encoder = createEncoder(encoderName, net, dim, rng.getInt());
	auto classifier = createClassifier(classifierName,
			encoder->getEdgeCodeDim(), rng.getInt());
	auto pr = std::make_shared<UECLPredictor<>>(net, encoder, classifier,
			rng.getInt());
	pr->setPosRatio(posRatio);
	pr->setNegRatio(negRatio);

	return predAll(pr);
}

void Predictor::predECL(std::vector<EdgeScore> &edgeScores,
		std::string encoderName, std::string classifierName, int dim,
		double posRatio, double negRatio, long int seed) {

	RandomGen rng(seed);
	auto encoder = createEncoder(encoderName, net, dim, rng.getInt());
	auto classifier = createClassifier(classifierName,
			encoder->getEdgeCodeDim(), rng.getInt());
	auto pr = std::make_shared<UECLPredictor<>>(net, encoder, classifier,
			rng.getInt());
	pr->setPosRatio(posRatio);
	pr->setNegRatio(negRatio);

	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopECL(int k, std::string encoderName,
		std::string classifierName, int dim, double posRatio, double negRatio,
		long int seed) {

	RandomGen rng(seed);
	auto encoder = createEncoder(encoderName, net, dim, rng.getInt());
	auto classifier = createClassifier(classifierName,
			encoder->getEdgeCodeDim(), rng.getInt());
	auto pr = std::make_shared<UECLPredictor<>>(net, encoder, classifier,
			rng.getInt());
	pr->setPosRatio(posRatio);
	pr->setNegRatio(negRatio);

	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllESM(std::string encoderName,
		std::string simMeasureName, int dim, long int seed) {

	RandomGen rng(seed);
	auto encoder = createEncoder(encoderName, net, dim, rng.getInt());
	auto simMeasure = createSimMeasure(simMeasureName);
	auto pr = std::make_shared<UESMPredictor<>>(net, encoder, simMeasure);

	return predAll(pr);
}

void Predictor::predESM(std::vector<EdgeScore> &edgeScores,
		std::string encoderName, std::string simMeasureName, int dim,
		long int seed) {

	RandomGen rng(seed);
	auto encoder = createEncoder(encoderName, net, dim, rng.getInt());
	auto simMeasure = createSimMeasure(simMeasureName);
	auto pr = std::make_shared<UESMPredictor<>>(net, encoder, simMeasure);

	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopESM(int k, std::string encoderName,
		std::string simMeasureName, int dim, long int seed) {

	RandomGen rng(seed);
	auto encoder = createEncoder(encoderName, net, dim, rng.getInt());
	auto simMeasure = createSimMeasure(simMeasureName);
	auto pr = std::make_shared<UESMPredictor<>>(net, encoder, simMeasure);

	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllFBM(int maxIter, long int seed) {
	auto pr = std::make_shared<UFBMPredictor<>>(net, seed);
	pr->setMaxIter(maxIter);
	return predAll(pr);
}

void Predictor::predFBM(std::vector<EdgeScore> &edgeScores, int maxIter,
		long int seed) {
	auto pr = std::make_shared<UFBMPredictor<>>(net, seed);
	pr->setMaxIter(maxIter);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopFBM(int k, int maxIter,
		long int seed) {
	auto pr = std::make_shared<UFBMPredictor<>>(net, seed);
	pr->setMaxIter(maxIter);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllHDI() {
	auto pr = std::make_shared<UHDIPredictor<>>(net);
	return predAll(pr);
}

void Predictor::predHDI(std::vector<EdgeScore> &edgeScores) {
	auto pr = std::make_shared<UHDIPredictor<>>(net);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopHDI(int k) {
	auto pr = std::make_shared<UHDIPredictor<>>(net);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllHPI() {
	auto pr = std::make_shared<UHPIPredictor<>>(net);
	return predAll(pr);
}

void Predictor::predHPI(std::vector<EdgeScore> &edgeScores) {
	auto pr = std::make_shared<UHPIPredictor<>>(net);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopHPI(int k) {
	auto pr = std::make_shared<UHPIPredictor<>>(net);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllHRG(int nbBeans, int nbSamples,
		long int seed) {
	auto pr = std::make_shared<UHRGPredictor<>>(net, seed);
	pr->setNbBeans(nbBeans);
	pr->setNbSamples(nbSamples);
	return predAll(pr);
}

void Predictor::predHRG(std::vector<EdgeScore> &edgeScores, int nbBeans,
		int nbSamples, long int seed) {
	auto pr = std::make_shared<UHRGPredictor<>>(net, seed);
	pr->setNbBeans(nbBeans);
	pr->setNbSamples(nbSamples);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopHRG(int k, int nbBeans, int nbSamples,
		long int seed) {
	auto pr = std::make_shared<UHRGPredictor<>>(net, seed);
	pr->setNbBeans(nbBeans);
	pr->setNbSamples(nbSamples);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllHYP(double m, double L, double gamma,
		double zeta, double T, long int seed) {
	auto pr = std::make_shared<UHYPPredictor<>>(net, seed);
	pr->setM(m);
	pr->setL(L);
	pr->setGamma(gamma);
	pr->setZeta(zeta);
	pr->setT(T);
	return predAll(pr);
}

void Predictor::predHYP(std::vector<EdgeScore> &edgeScores, double m, double L,
		double gamma, double zeta, double T, long int seed) {
	auto pr = std::make_shared<UHYPPredictor<>>(net, seed);
	pr->setM(m);
	pr->setL(L);
	pr->setGamma(gamma);
	pr->setZeta(zeta);
	pr->setT(T);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopHYP(int k, double m, double L,
		double gamma, double zeta, double T, long int seed) {
	auto pr = std::make_shared<UHYPPredictor<>>(net, seed);
	pr->setM(m);
	pr->setL(L);
	pr->setGamma(gamma);
	pr->setZeta(zeta);
	pr->setT(T);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllJID() {
	auto pr = std::make_shared<UJIDPredictor<>>(net);
	return predAll(pr);
}

void Predictor::predJID(std::vector<EdgeScore> &edgeScores) {
	auto pr = std::make_shared<UJIDPredictor<>>(net);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopJID(int k) {
	auto pr = std::make_shared<UJIDPredictor<>>(net);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllKAB(int horizLim) {
	auto pr = std::make_shared<UKABPredictor<>>(net);
	pr->setHorizLim(horizLim);
	return predAll(pr);
}

void Predictor::predKAB(std::vector<EdgeScore> &edgeScores, int horizLim) {
	auto pr = std::make_shared<UKABPredictor<>>(net);
	pr->setHorizLim(horizLim);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopKAB(int k, int horizLim) {
	auto pr = std::make_shared<UKABPredictor<>>(net);
	pr->setHorizLim(horizLim);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllLCP(double epsilon) {
	auto pr = std::make_shared<ULCPPredictor<>>(net);
	pr->setEpsilon(epsilon);
	return predAll(pr);
}

void Predictor::predLCP(std::vector<EdgeScore> &edgeScores, double epsilon) {
	auto pr = std::make_shared<ULCPPredictor<>>(net);
	pr->setEpsilon(epsilon);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopLCP(int k, double epsilon) {
	auto pr = std::make_shared<ULCPPredictor<>>(net);
	pr->setEpsilon(epsilon);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllLHN() {
	auto pr = std::make_shared<ULHNPredictor<>>(net);
	return predAll(pr);
}

void Predictor::predLHN(std::vector<EdgeScore> &edgeScores) {
	auto pr = std::make_shared<ULHNPredictor<>>(net);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopLHN(int k) {
	auto pr = std::make_shared<ULHNPredictor<>>(net);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllPAT() {
	auto pr = std::make_shared<UPATPredictor<>>(net);
	return predAll(pr);
}

void Predictor::predPAT(std::vector<EdgeScore> &edgeScores) {
	auto pr = std::make_shared<UPATPredictor<>>(net);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopPAT(int k) {
	auto pr = std::make_shared<UPATPredictor<>>(net);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllRAL() {
	auto pr = std::make_shared<URALPredictor<>>(net);
	return predAll(pr);
}

void Predictor::predRAL(std::vector<EdgeScore> &edgeScores) {
	auto pr = std::make_shared<URALPredictor<>>(net);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopRAL(int k) {
	auto pr = std::make_shared<URALPredictor<>>(net);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllRND(long int seed) {
	auto pr = std::make_shared<URNDPredictor<>>(net, seed);
	return predAll(pr);
}

void Predictor::predRND(std::vector<EdgeScore> &edgeScores, long int seed) {
	auto pr = std::make_shared<URNDPredictor<>>(net, seed);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopRND(int k, long int seed) {
	auto pr = std::make_shared<URNDPredictor<>>(net, seed);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllSAI() {
	auto pr = std::make_shared<USAIPredictor<>>(net);
	return predAll(pr);
}

void Predictor::predSAI(std::vector<EdgeScore> &edgeScores) {
	auto pr = std::make_shared<USAIPredictor<>>(net);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopSAI(int k) {
	auto pr = std::make_shared<USAIPredictor<>>(net);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllSBM(int maxIter, long int seed) {
	auto pr = std::make_shared<USBMPredictor<>>(net, seed);
	pr->setMaxIter(maxIter);
	return predAll(pr);
}

void Predictor::predSBM(std::vector<EdgeScore> &edgeScores, int maxIter,
		long int seed) {
	auto pr = std::make_shared<USBMPredictor<>>(net, seed);
	pr->setMaxIter(maxIter);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopSBM(int k, int maxIter,
		long int seed) {
	auto pr = std::make_shared<USBMPredictor<>>(net, seed);
	pr->setMaxIter(maxIter);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllSHP(long int seed) {
	auto pr = std::make_shared<USHPPredictor<>>(net, seed);
	return predAll(pr);
}

void Predictor::predSHP(std::vector<EdgeScore> &edgeScores, long int seed) {
	auto pr = std::make_shared<USHPPredictor<>>(net, seed);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopSHP(int k, long int seed) {
	auto pr = std::make_shared<USHPPredictor<>>(net, seed);
	return predTop(pr, k);
}

std::vector<EdgeScore> Predictor::predAllSOI() {
	auto pr = std::make_shared<USOIPredictor<>>(net);
	return predAll(pr);
}

void Predictor::predSOI(std::vector<EdgeScore> &edgeScores) {
	auto pr = std::make_shared<USOIPredictor<>>(net);
	pred(pr, edgeScores);
}

std::vector<EdgeScore> Predictor::predTopSOI(int k) {
	auto pr = std::make_shared<USOIPredictor<>>(net);
	return predTop(pr, k);
}

std::vector<EdgeScoreByID> Predictor::predAllADAByID() {
	auto pr = std::make_shared<UADAPredictor<>>(net);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllCNEByID() {
	auto pr = std::make_shared<UCNEPredictor<>>(net);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllCRAByID() {
	auto pr = std::make_shared<UCRAPredictor<>>(net);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllECLByID(std::string encoderName,
		std::string classifierName, int dim, double posRatio, double negRatio,
		long int seed) {

	RandomGen rng(seed);
	auto encoder = createEncoder(encoderName, net, dim, rng.getInt());
	auto classifier = createClassifier(classifierName,
			encoder->getEdgeCodeDim(), rng.getInt());
	auto pr = std::make_shared<UECLPredictor<>>(net, encoder, classifier,
			rng.getInt());
	pr->setPosRatio(posRatio);
	pr->setNegRatio(negRatio);

	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllESMByID(std::string encoderName,
		std::string simMeasureName, int dim, long int seed) {

	RandomGen rng(seed);
	auto encoder = createEncoder(encoderName, net, dim, rng.getInt());
	auto simMeasure = createSimMeasure(simMeasureName);
	auto pr = std::make_shared<UESMPredictor<>>(net, encoder, simMeasure);

	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllFBMByID(int maxIter,
		long int seed) {
	auto pr = std::make_shared<UFBMPredictor<>>(net, seed);
	pr->setMaxIter(maxIter);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllHDIByID() {
	auto pr = std::make_shared<UHDIPredictor<>>(net);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllHPIByID() {
	auto pr = std::make_shared<UHPIPredictor<>>(net);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllHRGByID(int nbBeans, int nbSamples,
		long int seed) {
	auto pr = std::make_shared<UHRGPredictor<>>(net, seed);
	pr->setNbBeans(nbBeans);
	pr->setNbSamples(nbSamples);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllHYPByID(double m, double L,
		double gamma, double zeta, double T, long int seed) {
	auto pr = std::make_shared<UHYPPredictor<>>(net, seed);
	pr->setM(m);
	pr->setL(L);
	pr->setGamma(gamma);
	pr->setZeta(zeta);
	pr->setT(T);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllJIDByID() {
	auto pr = std::make_shared<UJIDPredictor<>>(net);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllKABByID(int horizLim) {
	auto pr = std::make_shared<UKABPredictor<>>(net);
	pr->setHorizLim(horizLim);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllLCPByID(double epsilon) {
	auto pr = std::make_shared<ULCPPredictor<>>(net);
	pr->setEpsilon(epsilon);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllLHNByID() {
	auto pr = std::make_shared<ULHNPredictor<>>(net);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllPATByID() {
	auto pr = std::make_shared<UPATPredictor<>>(net);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllRALByID() {
	auto pr = std::make_shared<URALPredictor<>>(net);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllRNDByID(long int seed) {
	auto pr = std::make_shared<URNDPredictor<>>(net, seed);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllSAIByID() {
	auto pr = std::make_shared<USAIPredictor<>>(net);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllSBMByID(int maxIter,
		long int seed) {
	auto pr = std::make_shared<USBMPredictor<>>(net, seed);
	pr->setMaxIter(maxIter);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllSHPByID(long int seed) {
	auto pr = std::make_shared<USHPPredictor<>>(net, seed);
	return predAllByID(pr);
}

std::vector<EdgeScoreByID> Predictor::predAllSOIByID() {
	auto pr = std::make_shared<USOIPredictor<>>(net);
	return predAllByID(pr);
}

}
/* namespace Simp */
}
/* namespace LinkPred */

