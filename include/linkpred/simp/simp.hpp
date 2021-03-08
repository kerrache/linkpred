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
 * @brief Includes all headers of the essential interface.
 */

#ifndef SIMP_HPP_
#define SIMP_HPP_

/**
 * @defgroup Simp Simplified interface
 * This group contains a simplified interface for LinkPred that includes the essential functionalities.
 */

namespace LinkPred {
/**
 * @brief Simplified interface.
 * Contains a simplified interface for LinkPred that includes the essential functionalities.
 */
namespace Simp {
}
}

#include "edgescore.hpp"
#include "predictor.hpp"
#include "perfres.hpp"
#include "evaluator.hpp"

#endif /* SIMP_HPP_ */
