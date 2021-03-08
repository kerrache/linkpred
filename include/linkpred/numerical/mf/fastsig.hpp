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
 * @ingroup Numerical
 * @brief Contains the implementation of a fast sigmoid.
 */

#ifndef FASTSIG_HPP_
#define FASTSIG_HPP_

#include <type_traits>
#include <vector>
#include <cmath>

namespace LinkPred {

/**
 * @brief A fast sigmoid.
 * @tparam T Function value type (must be floating point).
 */
template<typename T> class FastSig {

	static_assert(std::is_floating_point<T>::value, "T must be floating point");

protected:
	int n = 1000; /**< Sigmoid table size. */
	T lowBound = -6; /**< Lower bound. */
	T uppBound = 6; /**< Upper bound. */
	std::vector<T> val; /**< Values. */

private:
	T interval; /**< uppBound - lowBound. */

	/**
	 * Initialize table.
	 */
	void initTab() {
		val.resize(n);
		interval = uppBound - lowBound;
		T x;
		for (int k = 0; k != n; k++) {
			x = interval * k / (n - 1) + lowBound;
			val[k] = 1.0 / (1 + std::exp(-x));
		}
	}

public:

	/**
	 * Default constructor.
	 */
	FastSig() {
		initTab();
	}

	/**
	 * Constructor.
	 * @param n Sigmoid table size.
	 * @param lowBound Lower bound. Any input smaller than lowBound is assigned the value 0.
	 * @param uppBound Upper bound. Any input greater than uppBound is assigned the value 1.
	 */
	FastSig(int n, T lowBound, T uppBound) :
			n(n), lowBound(lowBound), uppBound(uppBound) {
		initTab();
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	FastSig(FastSig const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	FastSig& operator =(FastSig const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	FastSig(FastSig &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	FastSig& operator =(FastSig &&that) = default;

	/**
	 * @param x Input value.
	 * @return Sigmoid(x) (approximated of course).
	 */
	T operator()(T x) {
		if (x > uppBound)
			return 1;
		else if (x < lowBound)
			return 0;
		int k = (x - lowBound) / interval * (n - 1);
		return val[k];
	}

	/**
	 * Destructor.
	 */
	virtual ~FastSig() = default;
};

} /* namespace LinkPred */

#endif /* FASTSIG_HPP_ */
