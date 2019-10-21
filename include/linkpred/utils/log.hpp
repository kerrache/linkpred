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
 * @brief Contains the implementation of a log class.
 */

#ifndef INCLUDE_LOG_HPP_
#define INCLUDE_LOG_HPP_

#include "linkpred/utils/loglevel.hpp"
#include <sstream>
#include <iostream>
#include <fstream>
#include <ctime>

namespace LinkPred {

/**
 * @brief A log class.
 */
class Log {

protected:
	std::ostringstream buffer; /**< The buffer to which the log is written. */

public:
	static constexpr std::ostream* DefaultOutputStream = &std::cerr; /**< Default output stream. */
	static const LogLevel DefaultLogLevel = logError; /**< Default log level. */
	static std::ostream *out; /**< Output file. */
	static LogLevel logLevel; /**< Log level. */

	/**
	 * Constructor.
	 */
	Log(LogLevel level = logError) {
		std::time_t t = std::time(NULL);
		char timeBuf[128];
		std::strftime(timeBuf, sizeof(timeBuf), "%F:%T", std::localtime(&t));
		buffer << level << " :"
				<< std::string(level > logDebug ? (level - logDebug) * 4 : 1,
						' ') << std::string(timeBuf) << " : ";
	}

	/**
	 * Copy constructor.
	 * @param value Value to be logged.
	 */
	template<typename T> Log & operator<<(T const & value) {
		buffer << value;
		return *this;
	}

	/**
	 * Destructor.
	 */
	~Log() {
		buffer << std::endl;
		(*out) << buffer.str();
	}
};

#ifndef NDEBUG
/**
 * Macro for writing to log.
 * @param level The message log level.
 * @param message The text that is added to the log.
 */
#define logger(level, message) \
if (level > Log::logLevel) ; \
else Log(level) << std::string(__FILE__)<< " in "<< std::string(__FUNCTION__) << ":" << __LINE__ << " : " << message;
#else
#define logger(level, message)
#endif

} /* namespace LinkPred */

#endif /* INCLUDE_LOG_HPP_ */
