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
 * @brief Contains the implementation of the dot product similarity measure.
 */

#ifndef DOTPROD_HPP_
#define DOTPROD_HPP_

#include "linkpred/ml/simmeasures/simmeasure.hpp"

namespace LinkPred {

/**
 * @brief A simple dot product similarity measure.
 */
class DotProd: public SimMeasure {

public:

	/**
	 * Constructor.
	 */
	DotProd() {
		name = "DTP";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	DotProd(DotProd const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	DotProd& operator =(DotProd const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	DotProd(DotProd &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	DotProd& operator =(DotProd &&that) = default;

	/**
	 * Compute the similarity between two vectors.
	 * @param v1 First vector.
	 * @param v2 Second vector. Must be of the same dimension as v1.
	 * @return The similarity between v1 and v2.
	 */
	virtual double sim(Vec const & v1, Vec const & v2);

	/**
	 * Destructor.
	 */
	virtual ~DotProd() = default;

};

} /* namespace LinkPred */

#endif /* DOTPROD_HPP_ */
