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

#include "linkpred/ml/simmeasures/pearson.hpp"
#include "linkpred/utils/log.hpp"
#include <stdexcept>
#include <cmath>

namespace LinkPred {


double Pearson::sim(Vec const & v1, Vec const & v2) {
	if (v1.size() != v2.size()) {
		throw std::invalid_argument(
				"The two vectors must have the same dimension");
	}

	double v1Bar = 0;
	for (int i = 0; i < v1.size(); i++) {
		v1Bar += v1[i];
	}
	v1Bar /= v1.size();

	double v2Bar = 0;
	for (int i = 0; i < v2.size(); i++) {
		v2Bar += v2[i];
	}
	v2Bar /= v2.size();

	double cov = 0;
	double var1 = 0;
	double var2 = 0;
	for (int i = 0; i < v1.size(); i++) {
		cov += (v1[i] - v1Bar) * (v2[i] - v2Bar);
		var1 += (v1[i] - v1Bar) * (v1[i] - v1Bar);
		var2 += (v2[i] - v2Bar) * (v2[i] - v2Bar);
	}

	return cov / (std::sqrt(var1 * var2));
}


} /* namespace LinkPred */
