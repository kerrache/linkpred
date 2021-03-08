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

#include "linkpred/ml/simmeasures/cosinesim.hpp"
#include "linkpred/utils/log.hpp"
#include <stdexcept>
#include <cmath>

namespace LinkPred {


double CosineSim::sim(Vec const & v1, Vec const & v2) {
	if (v1.size() != v2.size()) {
		throw std::invalid_argument(
				"The two vectors must have the same dimension");
	}

	double dp = 0;
	double nrmV1 = 0;
	double nrmV2 = 0;
	for (int i = 0; i < v1.size(); i++) {
		dp += v1[i] * v2[i];
		nrmV1 += v1[i] * v1[i];
		nrmV2 += v2[i] * v2[i];
	}

	return dp / (std::sqrt(nrmV1 * nrmV2));
}


} /* namespace LinkPred */
