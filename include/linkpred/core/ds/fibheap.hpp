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
 * @brief Contains the implementation of a templated Fibonacci heap.
 */

#ifndef CORE_FIBHEAP_HPP_
#define CORE_FIBHEAP_HPP_

#include "linkpred/utils/log.hpp"

#include <map>

namespace LinkPred {

/**
 * @brief A Fibonacci heap.
 * @tparam T The stored data type.
 * @tparam P The priority type.
 * @tparam ComparatorT The data comparator.
 * @tparam ComparatorP The priority comparator.
 */
template<typename T, typename P = int, typename ComparatorT = std::less<T>,
		typename ComparatorP = std::less<P>> class FibHeap {

protected:

	/**
	 * @brief A Fibonacci heap node.
	 * @tparam T The stored data type.
	 * @tparam P The priority type.
	 * @tparam ComparatorT The data comparator.
	 * @tparam ComparatorP The priority comparator.
	 */
	template<typename T, typename P = int, typename ComparatorT = std::less<T>,
			typename ComparatorP = std::less<P>> class Node {
	public:
		T data; /**< The node data. */
		P pr; /**< Priority. */
		std::size_t deg = 0; /**< The node degree. */
		bool marked = false; /**< Whether the node is marked. */
		Node<T, P, ComparatorT, ComparatorP>* parent = nullptr; /**< Pointer to the parent. */
		Node<T, P, ComparatorT, ComparatorP>* child = nullptr; /**< Pointer to the child. */
		Node<T, P, ComparatorT, ComparatorP>* left = this; /**< Pointer to the left node. */
		Node<T, P, ComparatorT, ComparatorP>* right = this; /**< Pointer to the right node. */

		/**
		 * Constructor.
		 * @param data The node data.
		 * @param pr The priority.
		 */
		Node(T const & data, P const & pr) :
				data(data), pr(pr) {
		}

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		Node(Node const & that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		Node & operator =(Node const & that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		Node(Node && that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		Node & operator =(Node && that) = default;

	};

	using NodeT = Node<T, P, ComparatorT, ComparatorP>; /**< Node type. */
	NodeT* minNode = nullptr; /**< Pointer to the kinimum node. */
	std::map<T, NodeT*, ComparatorT> idMap; /**< Map of elements to their nodes in the heap. */
	std::size_t n = 0; /**< Number of nodes. */

	/**
	 * @brief Merge two heaps.
	 * @param m1 The first heap.
	 * @param m2 The second heap.
	 * @return The merged heap.
	 */
	static NodeT* merge(NodeT* m1, NodeT* m2) {
		if (m1 == nullptr) {
			return m2;
		}
		if (m2 == nullptr) {
			return m1;
		}

		m2->left->right = m1->right;
		m1->right->left = m2->right;
		m1->right = m2;
		m2->left = m1;

		ComparatorP compare;
		if (compare(m2->pr, m1->pr)) {
			return m2;
		} else {
			return m1;
		}
	}

public:

	/**
	 * Constructor.
	 */
	FibHeap() = default;

	/**
	 * Push an element. If the element already exists, its priority is updated.
	 * @param elem The element to be pushed.
	 * @param pr The priority.
	 * @return True if the element is inserted, false if it already exists.
	 */
	bool push(T const & elem, P const & pr) {
		auto fit = idMap.find(elem);
		if (fit == idMap.end()) {
			auto p = new NodeT(elem, pr);
			if (minNode == nullptr) {
				minNode = p;
			} else {
				p->right = minNode->right;
				minNode->right = p;
				p->left = minNode;
				p->right->left = p;
				ComparatorP compare;

				if (compare(minNode->pr, p->pr)) {
					minNode = p;
				}
			}
			idMap[elem] = p;
			return true;
		} else {
			set(fit, pr);
			return false;
		}
	}

	/**
	 * Destructor.
	 */
	virtual ~FibHeap() = default;
};

} /* namespace LinkPred */

#endif /* INCLUDE_CORE_FIBHEAP_HPP_ */
