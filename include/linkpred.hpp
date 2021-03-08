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

/**
 * \file
 * @brief Includes all headers of the library.
 */

#ifndef LINKPRED_HPP_
#define LINKPRED_HPP_

/**
 * @mainpage LinkPred: A high performance library for link prediction in complex networks
 *
 * This is LinkPred, a high performance library for link prediction in complex networks.
 *
 * LinkPred provides the following functionalities:
 *
 *	Basic data structures to efficiently store and access network data.
 *
 *	Basic graph algorithms such graph traversal, shortest path algorithms, and graph embedding methods.
 *
 *	Implementation of several topological similarity index predictors, for example: common neighbors, Adamic-Adard index and Jackard index among other predictors (a full list is available in the library documentation).
 *
 *	Implementation of several state-of-the-art global link predictors (a full list is available in the library documentation).
 *
 *	Implementation of several link prediction algorithms based on graph embedding techniques.
 *
 *	Test data generation from ground truth networks.
 *
 *	Performance evaluation functionalities.
 */

/**
 * @brief Main namespace.
 * @details
 * The main namespace of the library.
 */
namespace LinkPred {
}

#include "LinkPredConfig.hpp"
#include "linkpred/core/core.hpp"
#include "linkpred/graphalg/graphalg.hpp"
#include "linkpred/predictors/predictors.hpp"
#include "linkpred/perf/perf.hpp"
#include "linkpred/utils/utils.hpp"
#include "linkpred/ml/ml.hpp"
#include "linkpred/numerical/numerical.hpp"
#include "linkpred/simp/simp.hpp"

#endif /* LINKPRED_HPP_ */
