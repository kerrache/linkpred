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
 * @brief Contains the implementation of a templated bidirectional half map.
 */

#ifndef BHMAP_HPP_
#define BHMAP_HPP_

#include <map>
#include <vector>

namespace LinkPred {

/**
 * @brief A bidirectional half map.
 * @details This is a special type of a bidirectional map, where user-provided keys are mapped to positions.
 * The latter are contiguous from 0 to n-1, and are assigned according to insertion order.
 * The lookup of position given the key is done in O(log n), whereas the reverse lookup is done is O(1).
 * @tparam K The key type.
 * @tparam P The position type (must be integer).
 * @tparam Comparator Key comparator.
 */
template<typename K, typename P = std::size_t,
		typename Comparator = std::less<K>> class Bhmap {
protected:
	std::map<K, P, Comparator> map; /**< Map mapping data to position. */
	std::vector<std::pair<P, K>> vec; /**< Vector of data. */

public:
	using k_const_iterator = typename std::map<K, P, Comparator>::const_iterator; /**< Constant iterator over keys. */
	using p_const_iterator = typename std::vector<std::pair<P, K>>::const_iterator; /**< Constant iterator over positions. */

	/**
	 * Constructor.
	 */
	Bhmap() = default;

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	Bhmap(Bhmap const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	Bhmap & operator =(Bhmap const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	Bhmap(Bhmap && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	Bhmap & operator =(Bhmap && that) = default;

	/**
	 * Add an element if it does not already exist.
	 * @param key The key.
	 * @return True if the key is added, false if it already exists.
	 */
	std::pair<P, bool> insert(K const & key) {
		auto fit = map.find(key);
		if (fit == map.end()) {
			P pos = static_cast<P>(map.size());
			map[key] = pos;
			vec.push_back(std::make_pair(pos, key));
			return std::make_pair(pos, true);
		} else {
			return std::make_pair(fit->second, false);
		}
	}

	/**
	 * @return The size of the map.
	 */
	std::size_t size() const {
		return map.size();
	}

	/**
	 * @return Begin iterator.
	 */
	p_const_iterator pbegin() const {
		return vec.begin();
	}

	/**
	 * @return End iterator.
	 */
	p_const_iterator pend() const {
		return vec.end();
	}

	/**
	 * @return Begin iterator.
	 */
	k_const_iterator kbegin() const {
		return map.begin();
	}

	/**
	 * @return End iterator.
	 */
	k_const_iterator kend() const {
		return map.end();
	}

	/**
	 * @return A reference to the mapped key.
	 */
	P const & pos(K const & key) const {
		return map.at(key);
	}

	/**
	 * @return A reference to the mapped key.
	 */
	K const & key(P const & pos) const {
		if ((pos >= 0) && (pos < static_cast<P>(vec.size()))) {
			return vec[pos].second;
		} else {
			throw std::out_of_range("Out of range position");
		}
	}

	/**
	 * Find a key. This operation is O(log(n)).
	 * @param key The key.
	 * @return An iterator to the key.
	 */
	k_const_iterator pfind(K const & key) const {
		return map.find(key);
	}

	/**
	 * Find a position. This operation is O(1).
	 * @param pos The position.
	 * @return An iterator to the position.
	 */
	p_const_iterator kfind(P const & pos) const {
		if ((pos >= 0) && (pos < static_cast<P>(vec.size()))) {
			return vec.begin() + pos;
		} else {
			return vec.end();
		}
	}

	/**
	 * Destructor.
	 */
	virtual ~Bhmap() = default;
};

} /* namespace LinkPred */

#endif /* BHMAP_HPP_ */
