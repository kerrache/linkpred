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

#include "linkpred/ml/simmeasures/lpsim.hpp"
#include "linkpred/utils/log.hpp"
#include <stdexcept>
#include <cmath>

namespace LinkPred {


double LPSim::sim(Vec const & v1, Vec const & v2) {
	if (v1.size() != v2.size()) {
		throw std::invalid_argument(
				"The two vectors must have the same dimension");
	}

	double lp = 0;
	for (int i = 0; i < v1.size(); i++) {
		lp += std::pow(v1[i] - v2[i], p);
	}
	return -std::pow(lp, 1.0 / p);
}


} /* namespace LinkPred */
