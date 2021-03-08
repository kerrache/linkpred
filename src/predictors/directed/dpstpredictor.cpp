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

#include "LinkPredConfig.hpp"
#include <linkpred/predictors/directed/dpstpredictor.hpp>
#include <linkpred/utils/miscutils.hpp>

namespace LinkPred {


template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void DPSTPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::loadEdgeScores(std::string fileName) {

	auto esv = Utils::template readEdgeScores<typename Network::Label>(
			fileName);
	edgeScores = net->template createEdgeMapSP<double>();
	for (auto it = esv.begin(); it != esv.end(); ++it) {
		auto i = it->i;
		auto j = it->j;
		double score = it->score;
		(*edgeScores)[net->makeEdge(net->getID(i), net->getID(j))] = score;
	}
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> double DPSTPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::score(Edge const & e) {
	return edgeScores->at(e);
}

#define DPSTPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef DPSTPREDICTOR_CPP


} /* namespace LinkPred */
