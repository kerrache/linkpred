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
 * @ingroup Core
 * @brief Contains the implementation of a templated map-priority queue with limit on the capacity.
 */

#ifndef LMAPQUEUE_HPP_
#define LMAPQUEUE_HPP_

#include "LinkPredConfig.hpp"
#include <map>
#include <queue>
#include <vector>
#include <iostream>
#include <cmath>

#ifdef LINKPRED_WITH_OPENMP
#include <omp.h>
#endif

namespace LinkPred {

/**
 * @brief A map-priority queue with limit on the capacity.
 * @details 
 * @tparam K The key type.
 * @tparam P The priority type.
 * @tparam KComparator Key comparator.
 * @tparam PComparator Priority comparator.
 */
template<typename K, typename P, typename KComparator = std::less<K>,
		typename PComparator = std::greater<P>> class LMapQueue {

protected:

	/**
	 * @brief Comparing pairs based on second element only.
	 */
	struct Comparator {
		/**
		 * Compares pairs based on second element.
		 * @param left First pair.
		 * @param right Second pair.
		 * @return True if left.second and righ.second comapre true.
		 */
		bool operator ()(std::pair<K, P> const & left,
				std::pair<K, P> const & right) {
			PComparator pcomp;
			return pcomp(left.second, right.second);
		}
	};
	std::map<K, P, KComparator> map; /**< Map mapping keys to priorities. */
	std::priority_queue<std::pair<K, P>, std::vector<std::pair<K, P>>,
			Comparator> pq; /**< Priority queue sorted according to priorities. */
	std::size_t l; /**< Limit capacity. */

public:

	/**
	 * Constructor.
	 * @param l The capacity limit of the container.
	 */
	LMapQueue(std::size_t l) :
			l(l) {

	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	LMapQueue(LMapQueue const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	LMapQueue & operator =(LMapQueue const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	LMapQueue(LMapQueue && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	LMapQueue & operator =(LMapQueue && that) = default;

	/**
	 * @return A reference to the element with the highest priority.
	 */
	std::pair<K, P> const & top() const {
		return pq.top();
	}

	/**
	 * Checks whether the container is empty.
	 * @return True if the container is empty, false otherwise.
	 */
	bool empty() const {
		return map.empty();
	}

	/**
	 * @return The number of elements in the container.
	 */
	std::size_t size() const {
		return map.size();
	}

	/**
	 * Print the content of the queue.
	 */
	void printPQ() {
		std::vector<std::pair<K, P>> vec;
		std::cout << "pq: ";
		while (!pq.empty()) {
			std::cout << "(" << pq.top().first << "," << pq.top().second
					<< ") ";
			vec.push_back(pq.top());
			pq.pop();
		}
		std::cout << std::endl;
		for (auto it = vec.begin(); it != vec.end(); ++it) {
			pq.push(*it);
		}
	}

	/**
	 * Insert an element into the map. Only the top l elements are kept.
	 * @param elem A pair where the first element is the key of the element top insert, and the second element is the priority of the element.
	 * @return True if the insertion takes place, false otherwise.
	 */
	bool push(std::pair<K, P> const & elem) {
		return push(elem.first, elem.second);
	}

	/**
	 * Insert an element into the map. Only the top l elements are kept.
	 * @param key The key of the element top insert.
	 * @param pr The priority of the element.
	 * @return True if the insertion takes place, false otherwise.
	 */
	bool push(K const & key, P const & pr) {
		//printPQ();
		PComparator cmp;
		auto fit = map.find(key);
		if (fit == map.end()) {
			if (pq.size() < l) {
				pq.push(std::make_pair(key, pr));
				map[key] = pr;
				return true;
			} else if (!cmp(pq.top().second, pr)) {
				map.erase(pq.top().first);
				pq.pop();
				pq.push(std::make_pair(key, pr));
				map[key] = pr;
				return true;
			}
		}
		return false;
	}

	/**
	 * Compares the input priority to the priority of the top element in the queue. The queue must not be empty.
	 * A true return value means that an element with priority pr will be inserted to the queue.
	 * @param pr The priority to be compared.
	 * @return True result of comparing the priority of top to pr.
	 */
	bool compareTop(P const & pr) {
		PComparator cmp;
		return !cmp(pq.top().second, pr);
	}

	/**
	 * 	Returns a reference to the priority of the element with key equivalent to key. If no such element exists, an exception of type std::out_of_range is thrown.
	 *  @param key The key of the element to find.
	 *  @return Reference to the priority of the requested element
	 */
	P & at(K const & key) {
		return map.at(key);
	}

	/**
	 * 	Returns a reference to the priority of the element with key equivalent to key. If no such element exists, an exception of type std::out_of_range is thrown.
	 *  @param key The key of the element to find.
	 *  @return Reference to the priority of the requested element
	 */
	P const & at(K const & key) const {
		return map.at(key);
	}

	/**
	 * Removes the top element from the priority queue.
	 */
	void pop() {
		map.erase(pq.top().first);
		pq.pop();
	}

	/**
	 * Removes all element with priority strictly less than the specified value.
	 * @param pr The threshold priority.
	 */
	void pop(P const pr) {
		PComparator cmp;
		while (!pq.empty()) {
			if (!cmp(pr, pq.top().second)) {
				break;
			}
			map.erase(pq.top().first);
			pq.pop();
		}
	}

	/**
	 * @return An iterator to the first element in the map.
	 */
	auto begin() {
		return map.begin();
	}

	/**
	 * @return A constant iterator to the first element in the map.
	 */
	auto begin() const {
		return map.begin();
	}

	/**
	 * @return A constant iterator to the first element in the map.
	 */
	auto cbegin() const {
		return map.cbegin();
	}

	/**
	 * @return An iterator to the one-past-the-last element in the map.
	 */
	auto end() {
		return map.end();
	}

	/**
	 * @return A constant iterator to the one-past-the-last element in the map.
	 */
	auto end() const {
		return map.end();
	}

	/**
	 * @return A constant iterator to the one-past-the-last element in the map.
	 */
	auto cend() const {
		return map.cend();
	}

	/**
	 * @return The value of the limiting capacity.
	 */
	std::size_t getL() const {
		return l;
	}

	/**
	 * Clear the content of the container.
	 */
	void clear() {
		map.clear();
		pq.clear();
	}

	/**
	 * @param key  key value of the elements to count.
	 * @return Number of elements with key that compares equivalent to key, which is either 1 or 0.
	 */
	std::size_t count(K const& key) const {
		return map.count(key);
	}

	/**
	 * @param key value of the element to search for.
	 * @return Iterator to an element with key equivalent to key. If no such element is found, past-the-end (see end()) iterator is returned.
	 */
	auto find(K const & key) const {
		return map.find(key);
	}

	/**
	 * Merge several queues by selecting at most the top k elements. The elements selected are unique. In case of duplication, the highest
	 * priority is considered. The input queues are emptied.
	 */
	template<typename RandomIterator> static LMapQueue<K, P, KComparator,
			PComparator> merge(std::size_t k, RandomIterator begin,
			RandomIterator end) {

		LMapQueue<K, P, KComparator, PComparator> mq(k);

		for (auto it = begin; it < end; ++it) {
			while (!it->empty()) {
				mq.push(it->top());
				it->pop();
			}
		}

		return mq;
	}

#ifdef LINKPRED_WITH_OPENMP
	// TODO Remove the condition that the range size must be equal to the number of threads
	/**
	 * Parallel merge of several queues by selecting at most the top k elements. The elements selected are unique. In case of duplication, the highest
	 * priority is considered. The input queues are emptied. The size of the range must be equal to the number of threads. Results are merged in the first
	 * queue. All other queues are destroyed in the process.
	 */
	template<typename RandomIterator> static void parMerge(RandomIterator begin,
			RandomIterator end) {
#pragma omp parallel
		{
			int nbThreads = omp_get_num_threads();
			int threadID = omp_get_thread_num();
			if (nbThreads != (end - begin)) {
				if (threadID == 0) {
					std::cerr << nbThreads << " " << (end - begin) << std::endl;
					throw std::invalid_argument(
							"The size of the range must be equal to the number of threads");
				}
			} else {
				for (int j = 0; j < std::ceil(std::log2(nbThreads)); j++) {
					int d = std::pow(2, j);
					int nd = std::pow(2, j + 1);
					if (threadID % nd == 0) { // Receiver
						if (threadID + d < nbThreads) {
							int src = threadID + d;
							while (!(begin + src)->empty()) {
								(begin + threadID)->push((begin + src)->top());
								(begin + src)->pop();
							}
						}
					}
					// Wait for all processes to finish
#pragma omp barrier
				}
			}
		}
	}
#endif

	/**
	 * Destructor.
	 */
	virtual ~LMapQueue() = default;
};

} /* namespace LinkPred */

#endif /* LMAPQUEUE_HPP_ */
