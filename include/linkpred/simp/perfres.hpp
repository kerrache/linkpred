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

/**
 * \file
 * @ingroup Simp
 * @brief Contains the definition of a structure to store performance results.
 */

#ifndef SIMPPERFRES_HPP_
#define SIMPPERFRES_HPP_

#include <string>

namespace LinkPred {
namespace Simp {

/**
 * @brief A structure to store performance results.
 */
struct PerfRes {
	std::string name; /**< Concatenation of the name of the performance mneasure and that of the predictor. */
	double res; /**< The result. */
};

}
/* namespace Simp */
}
/* namespace LinkPred */

#endif /* SIMPPERFRES_HPP_ */
