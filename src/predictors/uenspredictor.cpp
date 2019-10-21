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

#include <linkpred/predictors/uenspredictor.hpp>
#include "linkpred/utils/utilities.hpp"
#include "linkpred/utils/randomgen.hpp"
#include "linkpred/perf/perfevaluator.hpp"
#include <cmath>
#include <algorithm>

namespace LinkPred {

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UENSPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::init() {
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> double UENSPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::getPerf(
		typename std::vector<double>::iterator posBegin,
		typename std::vector<double>::iterator posEnd,
		typename std::vector<double>::iterator negBegin,
		typename std::vector<double>::iterator negEnd) const {

	SortOrder posSortOrder = SortOrder::None;
	SortOrder negSortOrder = SortOrder::None;
	PerfResults res;
	perfMeasure->eval(posBegin, posEnd, negBegin, negEnd, posSortOrder,
			negSortOrder, res);
	return res.at(perfMeasure->getName());
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UENSPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::selectEdgeScoreMethod() {

	logger(logDebug, "Selecting edge score method ...")

	std::vector<EdgeType> posLinks;
	posLinks.reserve(scoreMethodLinksRatio * net->getNbEdges());
	for (auto it = net->rndEdgesBegin(scoreMethodLinksRatio, rng.getInt());
			it != net->rndEdgesEnd(); ++it) {
		posLinks.push_back(*it);
	}
	std::vector<double> posLinksScore;
	posLinksScore.resize(posLinks.size());

	std::vector<EdgeType> negLinks;
	negLinks.reserve(scoreMethodLinksRatio * net->getNbNonEdges());
	for (auto it = net->rndNonEdgesBegin(scoreMethodLinksRatio, rng.getInt());
			it != net->rndNonEdgesEnd(); ++it) {
		negLinks.push_back(*it);
	}
	std::vector<double> negLinksScore;
	negLinksScore.resize(negLinks.size());

	if (perfMeasure == nullptr) {
		perfMeasure = std::make_shared<ROC<>>();
#ifdef WITH_OPENMP
		perfMeasure->setParallel(parallel);
#endif
	}

	EdgeScoreMethod bestMethod = AUT;
	double bestPerf = -1;

//	Utils::print(posLinks.begin(), posLinks.end(), "pos");
//	Utils::print(negLinks.begin(), negLinks.end(), "neg");

	{

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < negLinks.size(); i++) {
			negLinksScore[i] = ADAEdgeLength(negLinks[i]);
		}
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < posLinks.size(); i++) {
			posLinksScore[i] = ADAEdgeLength(posLinks[i]);
		}
//		Utils::print(posLinksScore.begin(), posLinksScore.end(), "pos");
//		Utils::print(negLinksScore.begin(), negLinksScore.end(), "neg");
		double perf = getPerf(posLinksScore.begin(), posLinksScore.end(),
				negLinksScore.begin(), negLinksScore.end());
		if (perf > bestPerf) {
			bestPerf = perf;
			bestMethod = ADA;
		}
		logger(logDebug1, "Perf of ADA : " << perf)
	}

	{

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < negLinks.size(); i++) {
			negLinksScore[i] = CNEEdgeLength(negLinks[i]);
		}
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < posLinks.size(); i++) {
			posLinksScore[i] = CNEEdgeLength(posLinks[i]);
		}
		double perf = getPerf(posLinksScore.begin(), posLinksScore.end(),
				negLinksScore.begin(), negLinksScore.end());
		if (perf > bestPerf) {
			bestPerf = perf;
			bestMethod = CNE;
		}
		logger(logDebug1, "Perf of CNE : " << perf)
	}

	{

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < negLinks.size(); i++) {
			negLinksScore[i] = HDIEdgeLength(negLinks[i]);
		}
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < posLinks.size(); i++) {
			posLinksScore[i] = HDIEdgeLength(posLinks[i]);
		}
		double perf = getPerf(posLinksScore.begin(), posLinksScore.end(),
				negLinksScore.begin(), negLinksScore.end());
		if (perf > bestPerf) {
			bestPerf = perf;
			bestMethod = HDI;
		}
		logger(logDebug1, "Perf of HDI : " << perf)
	}

	{

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < negLinks.size(); i++) {
			negLinksScore[i] = HPIEdgeLength(negLinks[i]);
		}
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < posLinks.size(); i++) {
			posLinksScore[i] = HPIEdgeLength(posLinks[i]);
		}
		double perf = getPerf(posLinksScore.begin(), posLinksScore.end(),
				negLinksScore.begin(), negLinksScore.end());
		if (perf > bestPerf) {
			bestPerf = perf;
			bestMethod = HPI;
		}
		logger(logDebug1, "Perf of HPI : " << perf)
	}

	{

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < negLinks.size(); i++) {
			negLinksScore[i] = JIDEdgeLength(negLinks[i]);
		}
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < posLinks.size(); i++) {
			posLinksScore[i] = JIDEdgeLength(posLinks[i]);
		}
		double perf = getPerf(posLinksScore.begin(), posLinksScore.end(),
				negLinksScore.begin(), negLinksScore.end());
		if (perf > bestPerf) {
			bestPerf = perf;
			bestMethod = JID;
		}
		logger(logDebug1, "Perf of JID : " << perf)
	}

	{

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < negLinks.size(); i++) {
			negLinksScore[i] = LHNEdgeLength(negLinks[i]);
		}
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < posLinks.size(); i++) {
			posLinksScore[i] = LHNEdgeLength(posLinks[i]);
		}
		double perf = getPerf(posLinksScore.begin(), posLinksScore.end(),
				negLinksScore.begin(), negLinksScore.end());
		if (perf > bestPerf) {
			bestPerf = perf;
			bestMethod = LHN;
		}
		logger(logDebug1, "Perf of LHN : " << perf)
	}

	{

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < negLinks.size(); i++) {
			negLinksScore[i] = PATEdgeLength(negLinks[i]);
		}
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < posLinks.size(); i++) {
			posLinksScore[i] = PATEdgeLength(posLinks[i]);
		}
		double perf = getPerf(posLinksScore.begin(), posLinksScore.end(),
				negLinksScore.begin(), negLinksScore.end());
		if (perf > bestPerf) {
			bestPerf = perf;
			bestMethod = PAT;
		}
		logger(logDebug1, "Perf of PAT : " << perf)
	}

	{

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < negLinks.size(); i++) {
			negLinksScore[i] = RALEdgeLength(negLinks[i]);
		}
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < posLinks.size(); i++) {
			posLinksScore[i] = RALEdgeLength(posLinks[i]);
		}
		double perf = getPerf(posLinksScore.begin(), posLinksScore.end(),
				negLinksScore.begin(), negLinksScore.end());
		if (perf > bestPerf) {
			bestPerf = perf;
			bestMethod = RAL;
		}
		logger(logDebug1, "Perf of RAL : " << perf)
	}

	{

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < negLinks.size(); i++) {
			negLinksScore[i] = SAIEdgeLength(negLinks[i]);
		}
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < posLinks.size(); i++) {
			posLinksScore[i] = SAIEdgeLength(posLinks[i]);
		}
		double perf = getPerf(posLinksScore.begin(), posLinksScore.end(),
				negLinksScore.begin(), negLinksScore.end());
		if (perf > bestPerf) {
			bestPerf = perf;
			bestMethod = SAI;
		}
		logger(logDebug1, "Perf of SAI : " << perf)
	}

	{

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < negLinks.size(); i++) {
			negLinksScore[i] = SOIEdgeLength(negLinks[i]);
		}
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < posLinks.size(); i++) {
			posLinksScore[i] = SOIEdgeLength(posLinks[i]);
		}
		double perf = getPerf(posLinksScore.begin(), posLinksScore.end(),
				negLinksScore.begin(), negLinksScore.end());
		if (perf > bestPerf) {
			bestPerf = perf;
			bestMethod = SOI;
		}
		logger(logDebug1, "Perf of SOI : " << perf)
	}

	if (bestMethod == AUT) {
		throw std::runtime_error("Unable to choose an edge length method");
	}

	scoreMethod = bestMethod;
//	std::cout << "Selected: " << EdgeScoreMethodToStr.at(scoreMethod)
//			<< std::endl;
	logger(logDebug, "Done")
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UENSPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::learn() {
	if (scoreMethod == AUT) {
		selectEdgeScoreMethod();
	}
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UENSPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::predict(EdgesRandomIteratorT begin,
		EdgesRandomIteratorT end, ScoresRandomIteratorT scores) {

	logger(logDebug, "Predicting links...")
	switch (scoreMethod) {

	case AUT: {
		throw std::logic_error(
				"Cannot predict with undetermined edge length method");
	}

	case CST: {
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = begin; it < end; ++it) {
			*(scores + (it - begin)) = CSTEdgeLength(*it);
		}
		break;
	}

	case ADA: {
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = begin; it < end; ++it) {
			*(scores + (it - begin)) = ADAEdgeLength(*it);
		}
		break;
	}

	case CNE: {
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = begin; it < end; ++it) {
			*(scores + (it - begin)) = CNEEdgeLength(*it);
		}
		break;
	}

	case HDI: {
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = begin; it < end; ++it) {
			*(scores + (it - begin)) = HDIEdgeLength(*it);
		}
		break;
	}

	case HPI: {
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = begin; it < end; ++it) {
			*(scores + (it - begin)) = HPIEdgeLength(*it);
		}
		break;
	}

	case JID: {
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = begin; it < end; ++it) {
			*(scores + (it - begin)) = JIDEdgeLength(*it);
		}
		break;
	}

	case LHN: {
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = begin; it < end; ++it) {
			*(scores + (it - begin)) = LHNEdgeLength(*it);
		}
		break;
	}

	case PAT: {
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = begin; it < end; ++it) {
			*(scores + (it - begin)) = PATEdgeLength(*it);
		}
		break;
	}

	case RAL: {
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = begin; it < end; ++it) {
			*(scores + (it - begin)) = RALEdgeLength(*it);
		}
		break;

	}

	case SAI: {
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = begin; it < end; ++it) {
			*(scores + (it - begin)) = SAIEdgeLength(*it);
		}
		break;
	}

	case SOI: {
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = begin; it < end; ++it) {
			*(scores + (it - begin)) = SOIEdgeLength(*it);
		}
		break;
	}

	default:
		throw std::runtime_error("Unknown edge length assignment strategy");
	}
	logger(logDebug, "Done.")
}

#define UENSPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef UENSPREDICTOR_CPP

}
/* namespace LinkPred */
