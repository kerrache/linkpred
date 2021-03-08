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
 * @brief Contains the definition of a structure to store the score of an edge.
 */

#ifndef SIMPEDGESCORE_HPP_
#define SIMPEDGESCORE_HPP_

#include <string>

namespace LinkPred {
namespace Simp {

/**
 * @brief A structure to store the score of an edge.
 */
struct EdgeScore {
	std::string i; /**< The label of the start node. */
	std::string j; /**< The label of the end node. */
	double score; /**< The score. */
};

/**
 * @brief A structure to store the score of an edge. The node IDs are used instead of labels.
 *
 * This is useful for reducing memory consumption for large networks.
 */
struct EdgeScoreByID {
	int i; /**< The ID of the start node. */
	int j; /**< The ID of the end node. */
	double score; /**< The score. */
};

}
/* namespace Simp */
}
/* namespace LinkPred */

#endif /* SIMPEDGESCORE_HPP_ */
