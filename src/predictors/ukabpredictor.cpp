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

#include "linkpred/predictors/ukabpredictor.hpp"
#include "linkpred/predictors/usnspredictor.hpp"
#include "linkpred/utils/utilities.hpp"
#include "linkpred/utils/randomgen.hpp"
#include "linkpred/perf/perfmeasure.hpp"
#include "linkpred/perf/networkmanipulator.hpp"
#include "linkpred/perf/perfevaluator.hpp"
#include <cmath>
#include <stdexcept>
#include <limits>
#ifdef WITH_OPENMP
#include <omp.h>
#endif
#ifdef WITH_MPI
#include <mpi.h>
#include <limits>
#include "linkpred/utils/utilities.hpp"
#endif

namespace LinkPred {

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UKABPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::init() {
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UKABPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::createLengthMap() {
	logger(logDebug, "Creating length map...")

	length = net->template createEdgeMapSP<double>();

	for (auto it = net->edgesBegin(); it < net->edgesEnd(); ++it) {
		(*length)[*it] = edgeLength(*it);
	}

	logger(logDebug, "Done")
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UKABPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::learn() {
	findLambda();
//	std::cout << "#lambda: " << lambda << std::endl;
	if (zeta == -1) {
		zeta = net->getCC();
		zeta = 1 / (1 + 0.5 * std::sqrt(1 - zeta) / std::sqrt(zeta));
	}
//	std::cout << "#zeta: " << zeta << std::endl;

	createLengthMap(); // Create the  edges length map
	distCalc = std::make_shared<ESPLDistCalculator<NetworkT>>(dijkstra, length,
			lim, cacheLevel);
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UKABPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::predict(EdgesRandomIteratorT begin,
		EdgesRandomIteratorT end, ScoresRandomIteratorT scores) {

	logger(logDebug, "Predicting links...")

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		*(scores + (it - begin)) = score(*it);
	}
	logger(logDebug, "Done")
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> double UKABPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::score(EdgeType const & e) {

	logger(logDebug, "Predicting links...")

	auto srcNode = NetworkT::start(e);
	auto endNode = NetworkT::end(e);
	auto dist = getDist(srcNode, endNode);
	return score(e, dist);
	logger(logDebug, "Done")
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> double UKABPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::score(EdgeType const & e, double dist) {

	logger(logDebug, "Predicting links...")
	auto i = NetworkT::start(e);
	auto j = NetworkT::end(e);

	double ki = net->getDeg(i);
	double kj = net->getDeg(j);
	switch (scoreT) {
	case Dist:
		return 1.0 / dist;
	case SumDist:
		return sum(ki, kj) / dist;
	case ProdDist:
		return prod(ki, kj) / dist;
	case SNSCSc:
		return ((1 - lambda) * sum(ki, kj) + lambda * nsc(i, j, ki, kj));
	case SNSCDist:
//		return std::pow(1 + sum(ki, kj), lambda)
//				* std::pow(1 + nsc(srcNode, endNode, ki, kj), 1 - lambda) / dist;

//		return sum(ki, kj) / dist;
//		std::cout << net->getLabel(i) << "\t" << net->getLabel(j) << "\t"
//				<< "sum: " << sum(ki, kj) << " nsc: " << nsc(i, j, ki, kj)
//				<< " dist: " << dist << std::endl;

		return ((1 - lambda) * sum(ki, kj) + lambda * nsc(i, j, ki, kj)) / dist;
//		return ((1 - lambda) * std::log2(1 + sum(ki, kj))
//				+ lambda * std::log2(1 + nsc(srcNode, endNode, ki, kj))) / dist;
//		return std::pow(sum(ki, kj), lambda)
//				* std::pow(nsc(srcNode, endNode, ki, kj), 1 - lambda) / dist;
	case RA1Dist:
		return RA1EdgeLength(e) / dist;

	case RA2Dist:
		return RA2EdgeLength(e) / dist;

	default:
		throw std::runtime_error("Unknown score type");
	}
	logger(logDebug, "Done")
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> std::size_t UKABPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::top(std::size_t k,
		EdgesRandomOutputIteratorT eit, ScoresRandomIteratorT sit) {

	std::vector<LMapQueue<EdgeType, double>> mqs;

#ifdef WITH_OPENMP
#pragma omp parallel if (parallel)
	{
#pragma omp critical(initLMapQueueArray)
		{
#endif
	mqs.push_back(LMapQueue<EdgeType, double>(k));
#ifdef WITH_OPENMP
		}
	}
#endif

// Computing loop start, end and step in case of distributed processing
#ifdef WITH_MPI
	int nbProcs;
	int procID;
	if (distributed) {
		MPI_Comm_size(comm, &nbProcs);
		MPI_Comm_rank(comm, &procID);
	} else {
		nbProcs = 1;
		procID = 0;
	}

	auto begin = net->nodesBegin() + (net->getNbNodes() / nbProcs) * procID
			+ std::min(static_cast<std::size_t>(procID),
					net->getNbNodes() % nbProcs);
	auto end = net->nodesBegin() + (net->getNbNodes() / nbProcs) * (procID + 1)
			+ std::min(static_cast<std::size_t>(procID + 1),
					net->getNbNodes() % nbProcs);
#else
	auto begin = net->nodesBegin();
	auto end = net->nodesEnd();
#endif

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel) schedule(dynamic) // We need a dynamic schedule because iteration vary greatly in duration
#endif
	for (auto nit = begin; nit < end; ++nit) {
#ifdef WITH_OPENMP
		int ind = omp_get_thread_num();
#else
		int ind = 0;
#endif

		auto i = nit->first;
		auto nnzSMap = getFiniteDistMap(i);
		for (auto it = nnzSMap->begin(); it != nnzSMap->end(); ++it) {
			auto j = it->first;
			if (i < j) {
				auto e = net->makeEdge(i, j);
				double sc = score(e, it->second.first);
				mqs[ind].push(e, sc);
			}
		}
	}

#ifdef WITH_OPENMP
		// Now we merge results obtained by all threads
		if (parallel) {
			LMapQueue<EdgeType, double>::parMerge(mqs.begin(), mqs.end());
		}
#endif

#ifdef WITH_MPI
	if (distributed) {
		// Now we merge over all processes
		Utilities::merge(mqs[0], comm);
	}
#endif

	for (auto it = mqs[0].begin(); it != mqs[0].end(); ++it) {
		*eit = it->first;
		++eit;
		*sit = it->second;
		++sit;
	}
	return mqs[0].size();
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UKABPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::findLambda() {

	switch (lambdaMeth) {
	case User:
		return;

	case Est: {
		double cc = net->getCC();
		lambda = std::log2(1 + std::sqrt(cc));
//		lambda = 1 / (1 + std::exp(-8 * lambda));
//		double n = net->getNbNodes();
//		double m = net->getNbEdges();
//		double de = 2 * m / (n * (n - 1));
//		de = std::log2(de);
//		de = 2 * de / (de - 1);
//		de = de * de / 4;
//		lambda *= de;
		return;
	}
	case Opt: {
		lambdaOpt();
		return;
	}
	case Scan: {
		lambdaScan();
		return;
	}
	default:
		throw std::runtime_error("Unknown lambda method");
	}
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UKABPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::lambdaOpt() {

	double r;

	std::vector<double> posSUM;
	std::vector<double> posNSC;
	r = std::min(1.0, (double) sampleSize / net->getNbEdges());
	for (auto it = net->rndEdgesBegin(r, rng.getInt());
			it != net->rndEdgesEnd(); ++it) {
		auto i = NetworkT::start(*it);
		auto j = NetworkT::end(*it);

		double ki = net->getDeg(i);
		double kj = net->getDeg(j);
		posSUM.push_back(sum(ki, kj));
		posNSC.push_back(nsc(i, j, ki, kj));
	}

	std::vector<double> negSUM;
	std::vector<double> negNSC;
	r = std::min(1.0, (double) sampleSize / net->getNbNonEdges());
	for (auto it = net->rndNonEdgesBegin(r, rng.getInt());
			it != net->rndNonEdgesEnd(); ++it) {
		auto i = NetworkT::start(*it);
		auto j = NetworkT::end(*it);

		double ki = net->getDeg(i);
		double kj = net->getDeg(j);
		negSUM.push_back(sum(ki, kj));
		negNSC.push_back(nsc(i, j, ki, kj));
	}

	CG::cg_stats stats;
	double lm;
	double tol = 1.0e-7;
//	std::size_t n = net->getNbNodes();
//	std::size_t m = net->getNbEdges();
//	double w = 1.0 - std::log2(1.0 + (2.0 * m) / (n * (n - 1.0)));
	double w = 0.5;
	KABLambdaCG pb(posSUM, negSUM, posNSC, negNSC, w);
	CG::CGDescent solver(&pb);
	solver.cg_descent(&lm, 1, &stats, nullptr, tol, nullptr, true);
	lambda = pb.getSol();

}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UKABPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::lambdaScan() {

	double posRatioParams = 0.1;
	double negRatioParams = 0.1;

#ifdef WITH_MPI
	int nbProcs;
	int procID;
	if (distributed) {
		MPI_Comm_size(comm, &nbProcs);
		MPI_Comm_rank(comm, &procID);
	} else {
		nbProcs = 1;
		procID = 0;
	}

	auto begin = lambdas.begin() + (lambdas.size() / nbProcs) * procID
	+ std::min(static_cast<std::size_t>(procID),
			lambdas.size() % nbProcs);
	auto end = lambdas.begin() + (lambdas.size() / nbProcs) * (procID + 1)
	+ std::min(static_cast<std::size_t>(procID + 1),
			lambdas.size() % nbProcs);

#else
	auto begin = lambdas.begin();
	auto end = lambdas.end();
#endif

	std::vector<double> cpm;

	for (auto it = begin; it < end; ++it) {
		cpm.push_back(0);
	}

	int nbTests = 1;
	for (int i = 0; i < nbTests; i++) {
		auto testData = NetworkManipulator<NetworkT>::createTestDataRem(net,
				posRatioParams, false, true, 0, false, negRatioParams,
				rng.getInt(), false);

		auto pm =
				std::make_shared<
						TPR<
								PredResults<
										TestData<NetworkT, std::vector<EdgeType>>,
										ULPredictor<NetworkT,
												typename std::vector<
														typename NetworkT::EdgeType>::const_iterator,
												ScoresRandomIteratorT,
												EdgesRandomOutputIteratorT>,
										std::vector<double>> >>(
						testData.getNbPos());
#ifdef WITH_OPENMP
		pm->setParallel(parallel);
#endif

		pm->setUseTopMethod(true);
		testData.genPos();
		//testData.genNeg(); // We do not need negatives, since we are using TPR
		testData.lock();

		PerfEvaluator<TestData<NetworkT, std::vector<EdgeType>>,
				ULPredictor<NetworkT,
						typename std::vector<typename NetworkT::EdgeType>::const_iterator,
						ScoresRandomIteratorT, EdgesRandomOutputIteratorT>,
				PredResults<TestData<NetworkT, std::vector<EdgeType>>,
						ULPredictor<NetworkT,
								typename std::vector<typename NetworkT::EdgeType>::const_iterator,
								ScoresRandomIteratorT,
								EdgesRandomOutputIteratorT>, std::vector<double>>,
				PerfMeasure<
						PredResults<TestData<NetworkT, std::vector<EdgeType>>,
								ULPredictor<NetworkT,
										typename std::vector<
												typename NetworkT::EdgeType>::const_iterator,
										ScoresRandomIteratorT,
										EdgesRandomOutputIteratorT>,
								std::vector<double>>>> perf(testData);

		perf.addPerfMeasure(pm);

		for (auto it = begin; it < end; ++it) {
			auto predictor =
					std::make_shared<
							USNSPredictor<NetworkT,
									typename std::vector<
											typename NetworkT::EdgeType>::const_iterator,
									ScoresRandomIteratorT,
									EdgesRandomOutputIteratorT>>(
							testData.getObsNet());
			predictor->setName(std::to_string(*it));
			predictor->setLambda(*it);
#ifdef WITH_OPENMP
			predictor->setParallel(parallel);
#endif
			perf.addPredictor(predictor);
		}

#ifdef WITH_OPENMP
		perf.setParallel(parallel);
#endif
		perf.eval();

		int k = 0;
		for (auto it = perf.resultsBegin(); it != perf.resultsEnd();
				++it, k++) {
			cpm[k] += it->second;
		}
	}

	// Compute average
	for (auto it = cpm.begin(); it != cpm.end(); ++it) {
		*it /= nbTests;
	}

	int bestLambdaInd = -1;
	double bestPerf = -1;
	for (auto it = cpm.begin(); it != cpm.end(); ++it) {
		double val = *it;
//		std::cout << "useSUM tpr: " << val << std::endl;
		if (val > bestPerf) {
			bestPerf = val;
			bestLambdaInd = static_cast<int>((it - cpm.begin())
					+ (begin - lambdas.begin()));
		}
	}

	// If distributed look for best value of w over all processors
#ifdef WITH_MPI
	if (distributed) {
		struct WValRank {
			double val;
			int rank;
		} in, out;
		in.val = bestPerf;
		in.rank = bestLambdaInd;
		MPI_Reduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, comm); // Processor 0 receives best i index
		bestLambdaInd = out.rank;
		bestPerf = out.val;
		MPI_Bcast(&bestLambdaInd, 1, MPI_INT, 0, comm); // Processor 0 broadcasts best index
		MPI_Barrier(comm);
	}
#endif
	lambda = lambdas[bestLambdaInd];
	std::cout << " # Estimated lambda :" << lambda << std::endl;
}

#define UKABPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef UKABPREDICTOR_CPP

}
/* namespace LinkPred */
