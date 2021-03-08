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
 * @ingroup ML
 * @brief Contains the interface of a similarity measure.
 */

#ifndef SIMMEASURE_HPP_
#define SIMMEASURE_HPP_

#include "linkpred/numerical/linear/vec.hpp"
#include <vector>

namespace LinkPred {

/**
 * @brief Interface of a similarity measure.
 */
class SimMeasure {

protected:
	std::string name; /**< The name of the similarity measure. */

public:

	/**
	 * Default constructor.
	 */
	SimMeasure() = default;

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	SimMeasure(SimMeasure const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	SimMeasure& operator =(SimMeasure const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	SimMeasure(SimMeasure &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	SimMeasure& operator =(SimMeasure &&that) = default;

	/**
	 * Compute the similarity between two vectors.
	 * @param v1 First vector.
	 * @param v2 Second vector. Must be of the same dimension as v1.
	 * @return The similarity between v1 and v2.
	 */
	virtual double sim(Vec const & v1,
			Vec const & v2) = 0;

	/**
	 * @return The name of the SimMeasure.
	 */
	const std::string& getName() const {
		return name;
	}

	/**
	 * Set the name of the SimMeasure.
	 * @param name The new name of the SimMeasure.
	 */
	void setName(const std::string &name) {
		this->name = name;
	}

	/**
	 * Destructor.
	 */
	virtual ~SimMeasure() = default;

};

} /* namespace LinkPred */

#endif /* SIMMEASURE_HPP_ */
