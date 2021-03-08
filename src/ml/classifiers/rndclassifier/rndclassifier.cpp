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

#include "LinkPredConfig.hpp"

#include "linkpred/ml/classifiers/rndclassifier/rndclassifier.hpp"
#include "linkpred/utils/log.hpp"

namespace LinkPred {


template<typename InRndIt, typename OutRndIt,
		typename ScoreRndIt> void RndClassifier<InRndIt,
		OutRndIt, ScoreRndIt>::learn(
		InRndIt trInBegin,
		InRndIt trInEnd,
		OutRndIt trOutBegin,
		OutRndIt trOutEnd) {

	logger(logDebug, "Learning...")

	logger(logDebug, "Done")
}

template<typename InRndIt, typename OutRndIt,
		typename ScoreRndIt> void RndClassifier<InRndIt,
		OutRndIt, ScoreRndIt>::predict(
		InRndIt testInputsBegin, InRndIt testInputsEnd,
		ScoreRndIt scoresBegin) {

	if (testInputsEnd - testInputsBegin == 0) {
		return;
	}

	auto sit = scoresBegin;
	for (auto it = testInputsBegin; it < testInputsEnd; ++it, ++sit) {
		*sit = rng.getDouble(0, 1);
	}
}

#define RNDCLASSIFIER_CPP
#include "linkpred/instantiations.hpp"
#undef RNDCLASSIFIER_CPP

} /* namespace LinkPred */

