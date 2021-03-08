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
 * @brief Contains the implementation of a random number generator.
 */

#ifndef RANDOMGEN_HPP_
#define RANDOMGEN_HPP_

#include <random>
#include <memory>

namespace LinkPred {

/**
 * @brief A random number generator.
 * @details This is mainly a wrapper that simplifies access to C++11 random generating classes/methods.
 */
class RandomGen {
	long int seed; /**< Seed. */
	std::mt19937_64 rng; /**< The random number generator (Mersenne-Twister). */

public:

	/**
	 * Constructor.
	 */
	RandomGen() {
		std::random_device rd;
		seed = rd();
		rng.seed(seed);
	}

	/**
	 * Constructor with seed.
	 * @param seed The seed.
	 */
	RandomGen(long int seed) :
			seed(seed), rng(seed) {
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	RandomGen(RandomGen const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	RandomGen & operator =(RandomGen const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	RandomGen(RandomGen && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	RandomGen & operator =(RandomGen && that) = default;

	/**
	 *@return The seed.
	 */
	auto getSeed() {
		return seed;
	}

	/**
	 * @return A random integer.
	 */
	auto getInt() {
		return rng();
	}

	/**
	 * @return A random unsigned integer in the interval [low, high] (inclusive of the boundaries).
	 * @param low The lower bound of the interval.
	 * @param high The upper bound of the interval.
	 */
	std::size_t getUInt(std::size_t low, std::size_t high) {
		std::uniform_int_distribution<std::size_t> dis(low, high);
		return dis(rng);
	}

	/**
	 * @return A random signed integer in the interval [low, high] (inclusive of the boundaries).
	 * @param low The lower bound of the interval.
	 * @param high The upper bound of the interval.
	 */
	int getSInt(int low, int high) {
		std::uniform_int_distribution<int> dis(low, high);
		return dis(rng);
	}

	/**
	 * Generate a uniformly distributed double in the specified interval.
	 * @param low The left limit of the interval.
	 * @param high The right limit of the interval.
	 * @return A double uniformly distributed in the interval [low, high).
	 */
	double getDouble(double low, double high) {
		std::uniform_real_distribution<double> dis(low, high);
		return dis(rng);
	}

	/**
	 * @return A uniformly distributed boolean (a fair coin toss).
	 */
	bool getBool() {
		std::uniform_int_distribution<int> dis(0, 1);
		return dis(rng) == 1;
	}

	/**
	 * @return Sample from a power law distribution.
	 * @param minV The minimum value.
	 * @param maxV The maximum value.
	 * @param gamma The exponent of the power law.
	 */
	double getPL(double minV, double maxV, double gamma) {
		// x = [(x1^(n+1) - x0^(n+1))*y + x0^(n+1)]^(1/(n+1))
		return std::pow(
				(std::pow(maxV, gamma + 1) - std::pow(minV, gamma + 1))
						* getDouble(0, 1) + std::pow(minV, gamma + 1),
				1 / (gamma + 1));
	}

	/**
	 * @param p The probability of success of the associated Bernouli distribution.
	 * @return A sample from a geometric distribution.
	 */
	unsigned long int getGeo(double p) {
		auto q = getDouble(0, 1);
		if (q < p) {
			return 0;
		} else {
			return static_cast<unsigned long int>(std::ceil(
					std::log(1 - q) / std::log(1 - p) - 1));
		}
	}

	/**
	 * Destructor.
	 */
	virtual ~RandomGen() = default;
};
} /* namespace LinkPred */

#endif /* RANDOMGEN_HPP_ */
