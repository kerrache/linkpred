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
 * @brief Contains the implementation of a templated binary heap.
 */

#ifndef BHEAP_HPP_
#define BHEAP_HPP_

#include "linkpred/utils/log.hpp"
#include <vector>
#include <map>
#include <iostream>

namespace LinkPred {

/**
 * @brief A binary heap.
 * @details This class implements a binary heap priority queue offering O(log n) insert, remove, and O(log n)^2 priority update.
 * @tparam T The stored data type.
 * @tparam P The priority type.
 * @tparam ComparatorT The data comparator.
 * @tparam ComparatorP The priority comparator.
 */
template<typename T, typename P = int, typename ComparatorT = std::less<T>,
		typename ComparatorP = std::less<P>> class BHeap {
protected:
	std::vector<std::pair<T, P>> heap; /**< The heap is stored as a vector. */
	std::map<T, std::size_t, ComparatorT> idMap; /**< Map of elements to their positions in the heap. */

	/**
	 * Heap down starting at a given index.
	 * @param i The index where the heap down procedure starts.
	 */
	void heapDown(std::size_t i) {
		ComparatorP compare;

		while (true) {
			std::size_t left = 2 * i;
			std::size_t right = 2 * i + 1;
			std::size_t selected = i;

			if ((left < heap.size())
					&& (compare(heap[left].second, heap[selected].second))) {
				selected = left;
			}

			if ((right < heap.size())
					&& (compare(heap[right].second, heap[selected].second))) {
				selected = right;
			}

			if (selected != i) {
				std::swap(heap[i], heap[selected]);
				idMap.at(heap[i].first) = i;
				idMap.at(heap[selected].first) = selected;
				i = selected;
			} else {
				break;
			}
		}
	}

	/**
	 * Heap up starting at a given index.
	 * @param i The index where the heap up procedure starts.
	 */
	void heapUp(std::size_t i) {
		ComparatorP compare;
		while (true) {
			std::size_t parent = i / 2;
			if ((parent > 0) && compare(heap[i].second, heap[parent].second)) {
				std::swap(heap[i], heap[parent]);
				idMap.at(heap[i].first) = i;
				idMap.at(heap[parent].first) = parent;
				i = parent;
			} else {
				break;
			}
		}
	}

	/**
	 * Set the priority of an element.
	 * @param fit Iterator pointing to the element.
	 * @param pr The new priority.
	 */
	void set(typename std::map<T, std::size_t, ComparatorT>::iterator fit,
			P const & pr) {
		auto i = fit->second;
		heap[i].second = pr;
		heapUp(i);
		heapDown(i);
	}

public:
	/**
	 * Constructor.
	 */
	BHeap() {
		heap.resize(1);
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	BHeap(BHeap const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	BHeap & operator =(BHeap const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	BHeap(BHeap && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	BHeap & operator =(BHeap && that) = default;

	/**
	 * @return The size of the heap.
	 */
	std::size_t size() {
		return heap.size() - 1;
	}

	/**
	 * Push an element. If the element already exists, its priority is updated.
	 * @param elem The element to be pushed.
	 * @param pr The priority.
	 * @return True if the element is inserted, false if it already exists.
	 */
	bool push(T const & elem, P const & pr) {
		logger(logDebug3, "Pushing " << elem << "\t" << pr << "...")
		auto fit = idMap.find(elem);
		if (fit == idMap.end()) {
			heap.push_back(std::make_pair(elem, pr));
			idMap[elem] = heap.size() - 1;
			heapUp(heap.size() - 1);
			logger(logDebug3, "Done")
			return true;
		} else {
			set(fit, pr);
			logger(logDebug3, "Done")
			return false;
		}
	}

	/**
	 * Return the element with the highest priority. This method does not remove the element.
	 * @return A constant reference to pair, where the first element is the data and the second is the associated priority.
	 */
	std::pair<T, P> const & top() const {
		return heap[1];
	}

	/**
	 * Remove the element with the highest priority.
	 */
	void pop() {
		logger(logDebug3,
				"Popping " + std::to_string(heap[1].first) + "\t"
						+ std::to_string(heap[1].second) + "...")
		idMap.erase(heap[1].first);
		if (heap.size() > 2) {
			heap[1] = std::move(heap[heap.size() - 1]);
			idMap.at(heap[1].first) = 1;
		}
		heap.pop_back();
		if (heap.size() > 2) {
			heapDown(1);
		}
		logger(logDebug3, "Done")
	}

	/**
	 * Set the priority of an element.
	 * @param elem The element.
	 * @param pr The new priority.
	 */
	void set(T const & elem, P const & pr) {
		auto i = idMap.at(elem);
		heap[i].second = pr;
		heapUp(i);
		heapDown(i);
	}

	/**
	 * Increase the priority of an element.
	 * @param elem The element.
	 * @param pr The new priority. This must not be smaller than the current one.
	 */
	void increase(T const & elem, P const & pr) {
		logger(logDebug3, "Heap increase " << elem << "\t" << pr)
		auto i = idMap.at(elem);
		heap[i].second = pr;
		heapUp(i);
	}

	/**
	 * Decrease the priority of an element.
	 * @param elem The element.
	 * @param pr The new priority. This must not be greater than the current one.
	 */
	void decrease(T const & elem, P const & pr) {
		auto i = idMap.at(elem);
		heap[i].second = pr;
		heapDown(i);
	}

	/**
	 * Try to increase the priority.
	 * @param elem The element.
	 * @param pr The new priority. The priority is updated only if this is greater than the current one.
	 * @return True if the priority is updated, false otherwise.
	 */
	bool tryIncrease(T const & elem, P const & pr) {
		ComparatorP compare;
		auto i = idMap.at(elem);
		if (compare(pr, heap[i].second)) {
			heap[i].second = pr;
			heapUp(i);
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Try to decrease the priority.
	 * @param elem The element.
	 * @param pr The new priority. The priority is updated only if this is smaller than the current one.
	 * @return True if the priority is updated, false otherwise.
	 */
	bool tryDecrease(T const & elem, P const & pr) {
		ComparatorP compare;
		auto i = idMap.at(elem);
		if (compare(heap[i].second, pr)) {
			heap[i].second = pr;
			heapDown(i);
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Check if an element exists.
	 * @param elem The element to be checked.
	 * @return True if \p elem exists in the heap, false otherwise.
	 */
	bool contains(T const & elem) {
		auto fit = idMap.find(elem);
		return fit != idMap.end();
	}

	/**
	 * Print the heap content to std::cout.
	 */
	void print() {
		std::cout << "Heap: " << std::endl;
		for (int i = 1; i < heap.size(); i++) {
			std::cout << "(" << heap[i].first << " , " << heap[i].second << ")"
					<< std::endl;
		}
		std::cout << "Map: " << std::endl;
		for (auto it = idMap.begin(); it != idMap.end(); ++it) {
			std::cout << it->first << " -> " << it->second << std::endl;
		}
		std::cout << std::endl;
	}

	/**
	 * Destructor.
	 */
	virtual ~BHeap() = default;
};

} /* namespace LinkPred */

#endif /* BHEAP_HPP_ */
