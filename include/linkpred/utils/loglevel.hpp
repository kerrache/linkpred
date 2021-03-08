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
 * @ingroup Utils
 * @brief Contains the definition of log levels.
 */

#ifndef LOGLEVEL_HPP_
#define LOGLEVEL_HPP_

namespace LinkPred {

/**
 * @brief Enumeration of log levels.
 */
enum LogLevel {
	logError, /**< At this level, only errors are reported. */
	logWarning, /**< Warnings are included at this level. */
	logInfo, /**< Running information included at this level. */
	logDebug, /**< Debug level. */
	logDebug1, /**< Debug level 1. */
	logDebug2, /**< Debug level 2. */
	logDebug3 /**< Debug level 3. */
};

} /* namespace LinkPred */

#endif /* INCLUDE_LOGLEVEL_HPP_ */

