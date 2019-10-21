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
 * @brief Contains the implementation of an undirected network data structure.
 */

#ifndef UNETWORK_HPP_
#define UNETWORK_HPP_

#include "linkpred/utils/log.hpp"
#include "linkpred/core/bhmap.hpp"
#include "linkpred/utils/randomgen.hpp"
#include <cmath>
#include <memory>
#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include <type_traits>
#include <string>

// TODO add iterators to edges and nonedges that iterate over all elements in a random order

namespace LinkPred {
/**
 * @brief This class represents an undirected network in the sense of graph theory.
 * @tparam LabelTypeT Type of external labels.
 * @tparam NodeIdTypeT Type of internal node IDs. This must be an unsigned integral type.
 * @tparam EdgeTypeT Type of edges. This must be an unsigned integral type having at least double the size of @p NodeIdTypeT \ref UNetwork.
 * @tparam Comparator Comparator of LabelTypeT used to store external labels in a map. 
 */
template<typename LabelTypeT = std::string, typename NodeIdTypeT = unsigned int,
		typename EdgeTypeT = unsigned long long int,
		typename Comparator = std::less<LabelTypeT>> class UNetwork {

	// Some static checks of the template types.
	static_assert(std::is_integral<NodeIdTypeT>::value, "NodeIdType must be integral");
	static_assert(std::is_unsigned<NodeIdTypeT>::value, "NodeIdType must be unsigned");
	static_assert(std::is_integral<EdgeTypeT>::value, "EdgeType must be integral");
	static_assert(std::is_unsigned<EdgeTypeT>::value, "EdgeType must be unsigned");
	static_assert(sizeof(EdgeTypeT) >= 2 * sizeof(NodeIdTypeT), "Condition sizeof(EdgeType) >= 2 * sizeof(NodeIdType) unsatisfied");

public:
	// Provided types.
	using LabelType = LabelTypeT; /**< External label type. */
	using NodeIdType = NodeIdTypeT; /**< Internal node ID type. */
	using EdgeType = EdgeTypeT; /**< Internal edge type. */
	using LabelIterator = typename Bhmap<LabelType, NodeIdType, Comparator>::k_const_iterator; /**< External node iterator that offers the mapping to internal IDs. */
	using NodeIterator = typename Bhmap<LabelType, NodeIdType, Comparator>::p_const_iterator; /**< Internal node iterator (random access iterator). */
	using EdgeIterator = typename std::vector<EdgeType>::const_iterator; /**< Edge iterator. */

protected:
	// Some constants used for manipulating edges.
	static constexpr unsigned int NbBitsShift = 8 * sizeof(NodeIdType); /**< Number of bits to be shifted to extract edge start node. */
	static constexpr EdgeType Mask = ((EdgeType) 1 << NbBitsShift) - 1; /**< Mask used to extract edge end node. */

public:

	/**
	 * @brief Node-degree iterator. This class can be used to iterate over pairs of node IDs and degrees.
	 */
	class NodeDegIterator: public std::iterator<std::random_access_iterator_tag,
			const std::pair<NodeIdType, std::size_t>, long int> {
		friend class UNetwork; /**< UNetwork is a friend. */

	protected:
		UNetwork const &net; /**< The network on which to iterate. */
		typename UNetwork::NodeIterator nit; /**< A node iterator. */
		std::pair<NodeIdType, std::size_t> nodeDeg; /**< To store current node degree. */

		/**
		 * Constructor.
		 * @param net The network on which to iterate.
		 * @param nit A node iterator.
		 */
		NodeDegIterator(UNetwork const &net,
				typename UNetwork::NodeIterator nit) :
				net(net), nit(nit) {
			nodeDeg = std::make_pair(0, 0);
		}

	public:
		using pointer = typename std::iterator<std::random_access_iterator_tag, const std::pair<NodeIdType, std::size_t>, long int>::pointer; /**< The pointer type associated with the iterator. */
		using reference = typename std::iterator<std::random_access_iterator_tag, const std::pair<NodeIdType, std::size_t>, long int>::reference; /**< The reference type associated with the iterator. */
		using difference_type = typename std::iterator<std::random_access_iterator_tag, const std::pair<NodeIdType, std::size_t>, long int>::difference_type; /**< The difference type associated with the iterator. */

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		NodeDegIterator(NodeDegIterator const &that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		NodeDegIterator& operator =(NodeDegIterator const &that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		NodeDegIterator(NodeDegIterator &&that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		NodeDegIterator& operator =(NodeDegIterator &&that) = default;

		/**
		 * Dereference operator.
		 * @return A reference to the object to which the iterator points.
		 */
		reference operator*() {
			nodeDeg = std::make_pair(nit->first, net.getDeg(nit->first));
			return nodeDeg;
		}

		/**
		 * Arrow operator.
		 * @return A pointer to the object to which the iterator points.
		 */
		pointer operator->() {
			nodeDeg = std::make_pair(nit->first, net.getDeg(nit->first));
			return &(nodeDeg);
		}

		/**
		 * Pre-increment operator.
		 * @return A reference to the new iterator.
		 */
		NodeDegIterator& operator++() {
			nit++;
			return *this;
		}

		/**
		 * Pre-decrement operator.
		 * @return A reference to the new iterator.
		 */
		NodeDegIterator& operator--() {
			nit--;
			return *this;
		}

		/**
		 * Post-increment operator.
		 * @return A reference to the new iterator.
		 */
		NodeDegIterator operator++(int) {
			auto that = *this;
			++(*this);
			return that;
		}

		/**
		 * Post-decrement operator.
		 * @return A reference to the new iterator.
		 */
		NodeDegIterator operator--(int) {
			auto that = *this;
			--(*this);
			return that;
		}

		/**
		 * Arithmetic + operator.
		 * @param n Increment value.
		 * @return The new iterator.
		 */
		NodeDegIterator operator+(const difference_type &n) const {
			return NodeDegIterator(net, nit + n);
		}

		/**
		 * Arithmetic += operator.
		 * @param n Increment value.
		 * @return A reference to the new iterator.
		 */
		NodeDegIterator& operator+=(const difference_type &n) {
			nit += n;
			return *this;
		}

		/**
		 * Arithmetic - operator.
		 * @param n Decrement value.
		 * @return The new iterator.
		 */
		NodeDegIterator operator-(const difference_type &n) const {
			return NodeDegIterator(net, nit - n);
		}

		/**
		 * Arithmetic -= operator.
		 * @param n Decrement value.
		 * @return A reference to the new iterator.
		 */
		NodeDegIterator& operator-=(const difference_type &n) {
			nit -= n;
			return *this;
		}

		/**
		 * Difference between the present iterator and the one passed as parameter.
		 * @param that The other iterator.
		 * @return The difference between the current and that iterator.
		 */
		difference_type operator-(const NodeDegIterator &that) const {
			return nit - that.nit;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this equals that.
		 */
		bool operator==(const NodeDegIterator &that) const {
			return nit == that.nit;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is not equal to that.
		 */
		bool operator!=(const NodeDegIterator &that) const {
			return !(*this == that);
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is less than that.
		 */
		bool operator<(const NodeDegIterator &that) const {
			return nit < that.nit;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is greater that.
		 */
		bool operator>(const NodeDegIterator &that) const {
			return nit > that.nit;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is greater or equal to that.
		 */
		bool operator<=(const NodeDegIterator &that) const {
			return !(*this > that);
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is less or equal to that.
		 */
		bool operator>=(const NodeDegIterator &that) const {
			return !(*this < that);
		}
	};

	/**
	 * @brief Nonedges iterator.
	 */
	class NonEdgeIterator: public std::iterator<std::random_access_iterator_tag,
			const EdgeType, long int> {

		friend class UNetwork; /**< UNetwork is a friend. */

	private:
		std::size_t nbNonEdges; /**< The number of negative edges. */

	protected:
		UNetwork const &net; /**< The network on which to iterate. */
		std::size_t ord; /**< The order of the current edge. */
		std::size_t j = 0; /**< A lower estimate of the index of first element in neja larger than ord. */
		bool jSynced; /**< Whether j is synced with ord. */
		EdgeType e = 0; /**< Non-edge pointed by the iterator. */

		/**
		 * Constructor.
		 * @param net The network on which to iterate.
		 * @param ord The order of the initial edge.
		 * @param j A lower estimate of the index of first element in neja larger than i.
		 */
		NonEdgeIterator(UNetwork const &net, std::size_t ord, std::size_t j) :
				net(net), ord(ord), j(j) {
			nbNonEdges = net.getNbNonEdges();
			jSynced = true;
		}

		/**
		 * Constructor.
		 * @param net The network on which to iterate.
		 * @param ord The order of the initial edge.
		 */
		NonEdgeIterator(UNetwork const &net, std::size_t ord) :
				net(net), ord(ord) {
			nbNonEdges = net.getNbNonEdges();
			jSynced = false;
		}

	public:
		using pointer = typename std::iterator<std::random_access_iterator_tag, const EdgeType, long int>::pointer; /**< The pointer type associated with the iterator. */
		using reference = typename std::iterator<std::random_access_iterator_tag, const EdgeType, long int>::reference; /**< The reference type associated with the iterator. */
		using difference_type = typename std::iterator<std::random_access_iterator_tag, const EdgeType, long int>::difference_type; /**< The difference type associated with the iterator. */

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		NonEdgeIterator(NonEdgeIterator const &that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		NonEdgeIterator& operator =(NonEdgeIterator const &that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		NonEdgeIterator(NonEdgeIterator &&that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		NonEdgeIterator& operator =(NonEdgeIterator &&that) = default;

		/**
		 * Dereference operator.
		 * @return A reference to the object to which the iterator points.
		 */
		reference operator*() {
			if (jSynced) {
				e = net.getNonEdge(ord, j);
			} else {
				e = net.getNonEdgeAndIndex(ord, j);
				jSynced = true;
			}
			return e;
		}

		/**
		 * Arrow operator.
		 * @return A pointer to the object to which the iterator points.
		 */
		pointer operator->() {
			if (jSynced) {
				e = net.getNonEdge(ord, j);
			} else {
				e = net.getNonEdgeAndIndex(ord, j);
				jSynced = true;
			}
			return &(e);
		}

		/**
		 * Pre-increment operator.
		 * @return A reference to the new iterator.
		 */
		NonEdgeIterator& operator++() {
			ord++;
			return *this;
		}

		/**
		 * Pre-decrement operator.
		 * @return A reference to the new iterator.
		 */
		NonEdgeIterator& operator--() {
			ord--;
			return *this;
		}

		/**
		 * Post-increment operator.
		 * @return A reference to the new iterator.
		 */
		NonEdgeIterator operator++(int) {
			auto that = *this;
			++(*this);
			return that;
		}

		/**
		 * Post-decrement operator.
		 * @return A reference to the new iterator.
		 */
		NonEdgeIterator operator--(int) {
			auto that = *this;
			--(*this);
			return that;
		}

		/**
		 * Arithmetic + operator.
		 * @param n Increment value.
		 * @return The new iterator.
		 */
		NonEdgeIterator operator+(const difference_type &n) const {
			return NonEdgeIterator(net, ord + n);
		}

		/**
		 * Arithmetic += operator.
		 * @param n Increment value.
		 * @return A reference to the new iterator.
		 */
		NonEdgeIterator& operator+=(const difference_type &n) {
			ord += n;
			jSynced = false;
			return *this;
		}

		/**
		 * Arithmetic - operator.
		 * @param n Decrement value.
		 * @return The new iterator.
		 */
		NonEdgeIterator operator-(const difference_type &n) const {
			return NonEdgeIterator(net, ord - n);
		}

		/**
		 * Arithmetic -= operator.
		 * @param n Decrement value.
		 * @return A reference to the new iterator.
		 */
		NonEdgeIterator& operator-=(const difference_type &n) {
			ord -= n;
			jSynced = false;
			return *this;
		}

		/**
		 * Difference between the present iterator and the one passed as parameter.
		 * @param that The other iterator.
		 * @return The difference between the current and that iterator.
		 */
		difference_type operator-(const NonEdgeIterator &that) const {
			return ord - that.ord;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this equals that.
		 */
		bool operator==(const NonEdgeIterator &that) const {
			return ord == that.ord;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is not equal to that.
		 */
		bool operator!=(const NonEdgeIterator &that) const {
			return !(*this == that);
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is less than that.
		 */
		bool operator<(const NonEdgeIterator &that) const {
			return ord < that.ord;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is greater that.
		 */
		bool operator>(const NonEdgeIterator &that) const {
			return ord > that.ord;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is greater or equal to that.
		 */
		bool operator<=(const NonEdgeIterator &that) const {
			return !(*this > that);
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is less or equal to that.
		 */
		bool operator>=(const NonEdgeIterator &that) const {
			return !(*this < that);
		}
	};

	/**
	 * @brief Randomized Nodes iterator.
	 * @details This is forward iterator that can be used to randomly sample
	 * a subset of the nodes.
	 */
	class RndNodeIterator: public std::iterator<std::input_iterator_tag,
			const std::pair<NodeIdType, LabelType>, long int> {
		friend class UNetwork;

	protected:
		typename UNetwork::NodeIterator nit; /**< Node iterator. */
		typename UNetwork::NodeIterator end; /**< The end of nodes. */
		double ratio; /**< Ratio of nodes to be selected. */
		RandomGen rng; /**< Random number generator. */

		/**
		 * Constructor.
		 * @param nit Starting node iterator.
		 * @param end The end of nodes.
		 * @param ratio The ratio of nodes to select.
		 * @param seed The random number generator's seed.
		 */
		RndNodeIterator(typename UNetwork::NodeIterator nit,
				typename UNetwork::NodeIterator end, double ratio,
				long int seed) :
				end(end), ratio(ratio), rng(seed) {
			this->nit = nit + rng.getGeo(ratio);
			if (this->nit > end) {
				this->nit = end;
			}
		}

		/**
		 * Constructor.
		 * @param nit Starting node iterator.
		 * @param end The end of nodes.
		 * @param ratio The ratio of nodes to select.
		 * @param rng A random number generator.
		 */
		RndNodeIterator(typename UNetwork::NodeIterator nit,
				typename UNetwork::NodeIterator end, double ratio,
				RandomGen rng) :
				end(end), ratio(ratio), rng(rng) {
			this->nit = nit + rng.getGeo(ratio);
			if (this->nit > end) {
				this->nit = end;
			}
		}

	public:
		using pointer = typename std::iterator<std::input_iterator_tag, const std::pair<NodeIdType, LabelType>, long int>::pointer; /**< The pointer type associated with the iterator. */
		using reference = typename std::iterator<std::input_iterator_tag, const std::pair<NodeIdType, LabelType>, long int>::reference;
		/**< The reference type associated with the iterator. */

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		RndNodeIterator(RndNodeIterator const &that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		RndNodeIterator& operator =(RndNodeIterator const &that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		RndNodeIterator(RndNodeIterator &&that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		RndNodeIterator& operator =(RndNodeIterator &&that) = default;

		/**
		 * Dereference operator.
		 * @return A reference to the object to which the iterator points.
		 */
		reference operator*() const {
			return *nit;
		}

		/**
		 * Arrow operator.
		 * @return A pointer to the object to which the iterator points.
		 */
		pointer operator->() const {
			return &(*nit);
		}

		/**
		 * Pre-increment operator.
		 * @return A reference to the new iterator.
		 */
		RndNodeIterator& operator++() {
			nit += (rng.getGeo(ratio) + 1);
			if (nit > end) {
				nit = end;
			}
			return *this;
		}

		/**
		 * Post-increment operator.
		 * @return A reference to the new iterator.
		 */
		RndNodeIterator operator++(int) {
			auto that = *this;
			++(*this);
			return that;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this equals that.
		 */
		bool operator==(const RndNodeIterator &that) const {
			return nit == that.nit;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is not equal to that.
		 */
		bool operator!=(const RndNodeIterator &that) const {
			return !(*this == that);
		}
	}
	;

	/**
	 * @brief Randomized edges iterator.
	 * @details This is forward iterator that can be used to randomly sample
	 * a subset of the edges.
	 */
	class RndEdgeIterator: public std::iterator<std::input_iterator_tag,
			const EdgeType, long int> {
		friend class UNetwork;

	protected:
		typename UNetwork::EdgeIterator eit; /**< Edge iterator. */
		typename UNetwork::EdgeIterator end; /**< The end of edges. */
		double ratio; /**< Ratio of edges to be selected. */
		RandomGen rng; /**< Random number generator. */

		/**
		 * Constructor.
		 * @param eit Starting edge iterator.
		 * @param end The end of edges.
		 * @param ratio The ratio of edges to select.
		 * @param seed The random number generator's seed.
		 */
		RndEdgeIterator(typename UNetwork::EdgeIterator eit,
				typename UNetwork::EdgeIterator end, double ratio,
				long int seed) :
				end(end), ratio(ratio), rng(seed) {
			this->eit = eit + rng.getGeo(ratio);
			if (this->eit > end) {
				this->eit = end;
			}
		}

		/**
		 * Constructor.
		 * @param eit Starting edge iterator.
		 * @param end The end of edges.
		 * @param ratio The ratio of edges to select.
		 * @param rng A random number generator.
		 */
		RndEdgeIterator(typename UNetwork::EdgeIterator eit,
				typename UNetwork::EdgeIterator end, double ratio,
				RandomGen rng) :
				end(end), ratio(ratio), rng(rng) {
			this->eit = eit + rng.getGeo(ratio);
			if (this->eit > end) {
				this->eit = end;
			}
		}

	public:
		using pointer = typename std::iterator<std::input_iterator_tag, const EdgeType, long int>::pointer; /**< The pointer type associated with the iterator. */
		using reference = typename std::iterator<std::input_iterator_tag, const EdgeType, long int>::reference; /**< The reference type associated with the iterator. */

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		RndEdgeIterator(RndEdgeIterator const &that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		RndEdgeIterator& operator =(RndEdgeIterator const &that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		RndEdgeIterator(RndEdgeIterator &&that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		RndEdgeIterator& operator =(RndEdgeIterator &&that) = default;

		/**
		 * Dereference operator.
		 * @return A reference to the object to which the iterator points.
		 */
		reference operator*() const {
			return *eit;
		}

		/**
		 * Arrow operator.
		 * @return A pointer to the object to which the iterator points.
		 */
		pointer operator->() const {
			return &(*eit);
		}

		/**
		 * Pre-increment operator.
		 * @return A reference to the new iterator.
		 */
		RndEdgeIterator& operator++() {
			eit += (rng.getGeo(ratio) + 1);
			if (eit > end) {
				eit = end;
			}
			return *this;
		}

		/**
		 * Post-increment operator.
		 * @return A reference to the new iterator.
		 */
		RndEdgeIterator operator++(int) {
			auto that = *this;
			++(*this);
			return that;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this equals that.
		 */
		bool operator==(const RndEdgeIterator &that) const {
			return eit == that.eit;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is not equal to that.
		 */
		bool operator!=(const RndEdgeIterator &that) const {
			return !(*this == that);
		}
	};

	/**
	 * @brief Randomized nonedges iterator.
	 * @details This is forward iterator that can be used to randomly sample
	 * a subset of the negative edges.
	 */
	class RndNonEdgeIterator: public std::iterator<std::input_iterator_tag,
			const EdgeType, long int> {
		friend class UNetwork;

	private:
		std::size_t nbNonEdges = 0; /**< NUmber of nonedges. */
	protected:
		UNetwork const &net; /**< Network on which to iterate. */
		EdgeType e = 0; /**< Current edge. */
		std::size_t ord; /**< Order of the current edge. */
		double ratio = 1; /**< Ratio of negative edges to select. */
		RandomGen rng; /**< A random number generator. */

		/**
		 * Constructor.
		 * @param net The network on which to iterate.
		 * @param ord The order of the initial edge.
		 */
		RndNonEdgeIterator(UNetwork const &net, std::size_t ord) :
				net(net), ord(ord) {
		}

		/**
		 * Constructor.
		 * @param net The network on which to iterate.
		 * @param ord The order of the initial edge.
		 * @param ratio The ratio of edges to select.
		 * @param seed The random number generator's seed.
		 */
		RndNonEdgeIterator(UNetwork const &net, std::size_t ord, double ratio,
				long int seed) :
				net(net), ratio(ratio), rng(seed) {
			nbNonEdges = net.getNbNonEdges();
			this->ord = ord;
			if (this->ord < nbNonEdges) {
				this->ord += rng.getGeo(ratio);
			}
			if (this->ord > nbNonEdges) {
				this->ord = nbNonEdges;
			}
		}

		/**
		 * Constructor.
		 * @param net The network on which to iterate.
		 * @param ord The order of the initial edge.
		 * @param ratio The ratio of edges to select.
		 * @param rng A random number generator.
		 */
		RndNonEdgeIterator(UNetwork const &net, std::size_t ord, double ratio,
				RandomGen rng) :
				net(net), ratio(ratio), rng(rng) {
			nbNonEdges = net.getNbNonEdges();
			this->ord = ord;
			if (this->ord < nbNonEdges) {
				this->ord += rng.getGeo(ratio);
			}
			if (this->ord > nbNonEdges) {
				this->ord = nbNonEdges;
			}
		}

	public:
		using pointer = typename std::iterator<std::input_iterator_tag, const EdgeType, long int>::pointer; /**< The pointer type associated with the iterator. */
		using reference = typename std::iterator<std::input_iterator_tag, const EdgeType, long int>::reference; /**< The reference type associated with the iterator. */

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		RndNonEdgeIterator(RndNonEdgeIterator const &that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		RndNonEdgeIterator& operator =(RndNonEdgeIterator const &that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		RndNonEdgeIterator(RndNonEdgeIterator &&that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		RndNonEdgeIterator& operator =(RndNonEdgeIterator &&that) = default;

		/**
		 * Dereference operator.
		 * @return A reference to the object to which the iterator points.
		 */
		reference operator*() {
			e = net.getNonEdge(ord);
			return e;
		}

		/**
		 * Arrow operator.
		 * @return A pointer to the object to which the iterator points.
		 */
		pointer operator->() {
			e = net.getNonEdge(ord);
			return &(e);
		}

		/**
		 * Pre-increment operator.
		 * @return A reference to the new iterator.
		 */
		RndNonEdgeIterator& operator++() {
			ord += (rng.getGeo(ratio) + 1);
			if (ord > nbNonEdges) {
				ord = nbNonEdges;
			}
			return *this;
		}

		/**
		 * Post-increment operator.
		 * @return A reference to the new iterator.
		 */
		RndNonEdgeIterator operator++(int) {
			auto that = *this;
			++(*this);
			return that;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this equals that.
		 */
		bool operator==(const RndNonEdgeIterator &that) const {
			return ord == that.ord;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is not equal to that.
		 */
		bool operator!=(const RndNonEdgeIterator &that) const {
			return !(*this == that);
		}

	};

	/**
	 * @brief A node map.
	 * @details This class can be used to assign a value to every node in the network.
	 * Access to values is done in constant time.
	 * @tparam ValueType Type of mapped values.
	 */
	template<typename ValueType> class NodeMap {
		friend class UNetwork;
	protected:
		UNetwork const &net; /**< The network. */
		std::vector<ValueType> values; /**< The values. */

		/**
		 * Constructor.
		 * @param net The network.
		 */
		NodeMap(UNetwork const &net) :
				net(net) {
			values.resize(net.getNbNodes());
		}

	public:

		using NodeMapIterator = typename std::vector<ValueType>::iterator; /**< Iterator on the map values. */
		using NodeMapConstIterator = typename std::vector<ValueType>::const_iterator; /**< A constant iterator on the map values. */

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		NodeMap(NodeMap const &that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		NodeMap& operator =(NodeMap const &that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		NodeMap(NodeMap &&that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		NodeMap& operator =(NodeMap &&that) = default;

		/**
		 * @param i A node ID.
		 * @return The value associated with the node i.
		 */
		inline ValueType operator[](NodeIdType const &i) const {
			return values[i];
		}

		/**
		 * @param i A node ID.
		 * @return A reference to the value associated with the node i.
		 */
		inline ValueType& operator[](NodeIdType const &i) {
			return values[i];
		}

		/**
		 * @param i A node ID.
		 * @return The value associated with the node i.
		 */
		inline ValueType at(NodeIdType const &i) const {
			return values.at(i);
		}

		/**
		 * @return An iterator to the first element in the map.
		 */
		NodeMapIterator begin() {
			return values.begin();
		}

		/**
		 * @return An iterator to one-past-the-last element in the map.
		 */
		NodeMapIterator end() {
			return values.end();
		}

		/**
		 * @return A constant iterator to the first element in the map.
		 */
		NodeMapConstIterator cbegin() const {
			return values.cbegin();
		}

		/**
		 * @return A constant iterator to one-past-the-last element in the map.
		 */
		NodeMapConstIterator cend() const {
			return values.cend();
		}

		/**
		 * Destructor.
		 */
		~NodeMap() = default;
	};

	template<typename ValueT> using NodeMapSP = std::shared_ptr<NodeMap<ValueT>>;
	/**< Shared pointer to a node map. */

	/**
	 * @brief A sparse node map.
	 * @details This class can be used to assign a value to every node in the network.
	 * Access to values is done in constant time.
	 * @tparam ValueType Type of mapped values.
	 */
	template<typename ValueType> class NodeSMap {
		friend class UNetwork;
	protected:
		UNetwork const &net; /**< The network. */
		std::map<NodeIdType, ValueType> values; /**< The values. */
		ValueType defVal; /**< Default value. */

		/**
		 * Constructor.
		 * @param net The network.
		 * @param defVal Default value.
		 */
		NodeSMap(UNetwork const &net, ValueType const defVal) :
				net(net), defVal(defVal) {
		}

	public:

		using NodeSMapIterator = typename std::map<NodeIdType, ValueType>::iterator; /**< Iterator on the map values. */
		using NodeSMapConstIterator = typename std::map<NodeIdType, ValueType>::const_iterator; /**< A constant iterator on the map values. */

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		NodeSMap(NodeSMap const &that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		NodeSMap& operator =(NodeSMap const &that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		NodeSMap(NodeSMap &&that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		NodeSMap& operator =(NodeSMap &&that) = default;

		/**
		 * @param i A node ID.
		 * @return An iterator to the element with key i.
		 */
		auto find(NodeIdType const &i) {
			return values.find(i);
		}

		/**
		 * @param i A node ID.
		 * @return The value associated with the node i.
		 */
		inline ValueType operator[](NodeIdType const &i) const {
			auto fit = values.find(i);
			if (fit == values.end()) {
				return defVal;
			} else {
				return fit->second;
			}
		}

		/**
		 * @param i A node ID.
		 * @return A reference to the value associated with the node i.
		 */
		inline ValueType& operator[](NodeIdType const &i) {
			return values[i];
		}

		/**
		 * @param i A node ID.
		 * @return The value associated with the node i.
		 */
		inline ValueType at(NodeIdType const &i) const {
			auto fit = values.find(i);
			if (fit == values.end()) {
				return defVal;
			} else {
				return fit->second;
			}
		}

		/**
		 * @return An iterator to the first element in the map.
		 */
		NodeSMapIterator begin() {
			return values.begin();
		}

		/**
		 * @return An iterator to one-past-the-last element in the map.
		 */
		NodeSMapIterator end() {
			return values.end();
		}

		/**
		 * @return A constant iterator to the first element in the map.
		 */
		NodeSMapConstIterator cbegin() const {
			return values.cbegin();
		}

		/**
		 * @return A constant iterator to one-past-the-last element in the map.
		 */
		NodeSMapConstIterator cend() const {
			return values.cend();
		}

		/**
		 * @return The number of elements in the map (those that were explicitly inserted).
		 */
		std::size_t size() const {
			return values.size();
		}

		/**
		 * Destructor.
		 */
		~NodeSMap() = default;
	};

	template<typename ValueT> using NodeSMapSP = std::shared_ptr<NodeSMap<ValueT>>;
	/**< Shared pointer to a node map. */

	/**
	 * @brief An edge map.
	 * @details This class can be used to assign a value to every edge in the network.
	 * Access to values is done in logarithmic time.
	 * @tparam ValueType Type of mapped values.
	 */
	template<typename ValueType> class EdgeMap {
		friend class UNetwork;
	protected:
		UNetwork const &net; /**< The network. */
		std::map<EdgeType, ValueType> values; /**< The values. */

		/**
		 * Constructor.
		 * @param net The network.
		 */
		EdgeMap(UNetwork const &net) :
				net(net) {
		}

	public:

		using EdgeMapIterator = typename std::map<EdgeType, ValueType>::iterator; /**< Iterator on the map values. */
		using EdgeMapConstIterator = typename std::map<EdgeType,ValueType>::const_iterator; /**< A constant iterator on the map values. */

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		EdgeMap(EdgeMap const &that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		EdgeMap& operator =(EdgeMap const &that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		EdgeMap(EdgeMap &&that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		EdgeMap& operator =(EdgeMap &&that) = default;

		/**
		 * @param e An edge.
		 * @return The value associated with the edge e.
		 */
		inline ValueType operator[](EdgeType const &e) const {
			auto fit = values.find(e);
			if (fit != values.end()) {
				return fit->second;
			}
			return values.at(UNetwork::reverseEdge(e));
		}

		/**
		 * @param e An edge.
		 * @return A reference to the value associated with the edge e.
		 */
		inline ValueType& operator[](EdgeType const &e) {
			auto fit = values.find(e);
			if (fit != values.end()) {
				return fit->second;
			}
			return values[UNetwork::reverseEdge(e)];
		}

		/**
		 * @param e An edge.
		 * @return The value associated with the edge e.
		 */
		inline ValueType at(EdgeType const &e) const {
			auto fit = values.find(e);
			if (fit != values.end()) {
				return fit->second;
			}
			return values.at(UNetwork::reverseEdge(e));
		}

		/**
		 * @return An iterator to the first element in the map.
		 */
		EdgeMapIterator begin() {
			return values.begin();
		}

		/**
		 * @return An iterator to one-past-the-last element in the map.
		 */
		EdgeMapIterator end() {
			return values.end();
		}

		/**
		 * @return A constant iterator to the first element in the map.
		 */
		EdgeMapConstIterator cbegin() const {
			return values.cbegin();
		}

		/**
		 * @return A constant iterator to one-past-the-last element in the map.
		 */
		EdgeMapConstIterator cend() const {
			return values.cend();
		}

		/**
		 * @return The size of the map.
		 */
		std::size_t size() {
			return values.size();
		}

		/**
		 * Destructor.
		 */
		~EdgeMap() = default;
	};

	template<typename ValueT> using EdgeMapSP = std::shared_ptr<EdgeMap<ValueT>>;
	/**< Shared pointer to an edge map. */

	/**
	 * Make an edge in internal representation out of two nodes' internal IDs.
	 * @param i The starting node.
	 * @param j The end node.
	 * @return The edge (i, j).
	 */
	inline static EdgeType makeEdge(NodeIdType const &i, NodeIdType const &j) {
		return ((EdgeType) i << NbBitsShift) + j;
	}

	/**
	 * @param e An edge (i,j).
	 * @return The edge (j, i).
	 */
	inline static EdgeType reverseEdge(EdgeType const &e) {
		auto i = start(e);
		auto j = end(e);
		return ((EdgeType) j << NbBitsShift) + i;
	}

	/**
	 *@param edge An edge.
	 *@return The starting node of edge.
	 */
	inline static NodeIdType const start(EdgeType const &edge) {
		return static_cast<NodeIdType>(edge >> NbBitsShift);
	}

	/**
	 *@param edge An edge.
	 *@return The end node of edge.
	 */
	inline static NodeIdType const end(EdgeType const &edge) {
		return static_cast<NodeIdType>(edge & Mask);
	}

	/**
	 * Compare edge ends.
	 * @param e1 First edge.
	 * @param e2 Second edge.
	 * @return True if the end of e1 is smaller than that of e2 (comparison is based on node IDs).
	 */
	inline static bool compareEdgeEnd(EdgeType const &e1, EdgeType const &e2) {
		return (e1 & Mask) < (e2 & Mask);
	}

protected:
	//protected:
	std::size_t nbNodes = 0; /**< Number of nodes. */
	std::size_t nbEdges = 0; /**< Number of edges. */
	std::size_t nbNonEdges = 0; /**< Number of non-edges. */
	Bhmap<LabelType, NodeIdType> nodesIdMap; /**< A bhmap mapping external node IDs with internal indexes. */
	std::vector<std::size_t> sia; /**< Adjacency matrix row index. */
	std::vector<EdgeType> sja; /**< Adjacency matrix column index. */
	std::vector<std::size_t> aia; /**< Adjacency matrix row index. */
	std::vector<EdgeType> aja; /**< Edges where i < j (no repetition of edges (i,j) and (j,i) ) . */
	std::vector<EdgeType> neja; /**< Used to iterate over negative edges. */
	bool assembled = false; /**< To indicate whether the network has been assembled or not. */
	std::size_t minDeg = 0; /**< Minimum degree. */
	std::size_t maxDeg = 0; /**< Maximum degree. */
	double avgDeg = 0; /**< Average degree. */

	/**
	 * @param i Order of the non-edge.
	 * @param j An estimate of the index of first element in neja larger than i.
	 * @return A non-edge of order i.
	 */
	inline EdgeType getNonEdge(std::size_t const &i, std::size_t &j) const {
		// The following two loops give the same result as std::upper_bound(neja.begin(), neja.end(), i)
		while (neja[j] > i && j > 0) {
			j--;
		}
		while (neja[j] <= i) {
			j++;
		}
		if (j == 0) {
			return coupleAtOrd(i);
		} else {
			return coupleAtOrd(coupleOrd(aja[j - 1]) + i - neja[j - 1] + 1);
		}
	}

	/**
	 * @param i Order of the non-edge.
	 * @return A non-edge of order i.
	 */
	inline EdgeType getNonEdge(std::size_t const &i) const {
		auto uit = std::upper_bound(neja.begin(), neja.end(), i);
		auto j = uit - neja.begin();
		if (j == 0) {
			return coupleAtOrd(i);
		} else {
			return coupleAtOrd(coupleOrd(aja[j - 1]) + i - neja[j - 1] + 1);
		}
	}

	/**
	 * @param i Order of the non-edge.
	 * @param j As output only. It is set to the index of first element in neja larger than i.
	 * @return A non-edge of order i.
	 */
	inline EdgeType getNonEdgeAndIndex(std::size_t const &i,
			std::size_t &j) const {
		auto uit = std::upper_bound(neja.begin(), neja.end(), i);
		j = static_cast<std::size_t>(uit - neja.begin());
		if (j == 0) {
			return coupleAtOrd(i);
		} else {
			return coupleAtOrd(coupleOrd(aja[j - 1]) + i - neja[j - 1] + 1);
		}
	}

	/**
	 * Increments an edge.
	 * @param e The edge to be incremented.
	 * @return The next edge.
	 */
	inline EdgeType increment(EdgeType const &e) const {
		auto st = start(e);
		auto en = end(e);
		en++;
		if (en == nbNodes) {
			st++;
			en = st + 1;
		}
		return makeEdge(st, en);
	}

	/**
	 * Increments an edge.
	 * @param e The edge to be decremented.
	 * @return The previous edge.
	 */
	EdgeType decrement(EdgeType const &e) const {
		auto st = start(e);
		auto en = end(e);
		if (en == st + 1) {
			st--;
			en = nbNodes - 1;
		} else {
			en--;
		}
		return makeEdge(st, en);
	}

private:
	std::set<EdgeType> edgeSet; /**< Temporary edge list used before assembly. */

public:

	/**
	 * Default constructor.
	 */
	UNetwork() = default;

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	UNetwork(UNetwork const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	UNetwork& operator =(UNetwork const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	UNetwork(UNetwork &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	UNetwork& operator =(UNetwork &&that) = default;

	/**
	 * Add a node.
	 * @param nodeId The ID of the node.
	 * @return An std::pair, where first is the internal ID, and second is a boolean which is true if the node is actually added.
	 */
	std::pair<NodeIdType, bool> addNode(LabelType const &nodeId);

	/**
	 * Translates from external label to internal IDs. This method is O(log n), where n is the number of nodes.
	 * @param label An external node label.
	 * @return The internal ID of label;
	 */
	inline NodeIdType getID(LabelType const &label) const {
		return nodesIdMap.pos(label);
	}

	/**
	 * @param iid An internal node ID.
	 * @return The external label of the node iid.
	 */
	inline LabelType getLabel(NodeIdType const &iid) const {
		return nodesIdMap.key(iid);
	}

	/**
	 * @param label An external node ID.
	 * @return Iterator to the external node ID.
	 */
	LabelIterator findLabel(LabelType const &label) const {
		return nodesIdMap.pfind(label);
	}

	/**
	 * @param iid An internal node ID.
	 * @return Iterator to the internal node ID.
	 */
	NodeIterator findNode(NodeIdType const &iid) const {
		return nodesIdMap.kfind(iid);
	}

	/**
	 * Add an edge.
	 * @param i The starting node.
	 * @param j The end node.
	 */
	void addEdge(NodeIdType const &i, NodeIdType const &j);

	/**
	 *@param iid An internal node ID.
	 *@return An iterator to the first neighbor of iid.
	 */
	EdgeIterator neighborsBegin(NodeIdType const &iid) const {
		return sja.begin() + sia[iid];
	}

	/**
	 * @param e An edge.
	 * @return The order of the edge
	 */
	inline std::size_t coupleOrd(EdgeType const &e) const {
		std::size_t i = start(e);
		std::size_t j = end(e);
		return i * (2 * nbNodes - i - 3) / 2 + j - 1;
	}

	/**
	 * @param ord The order of an edge.
	 * @return The edge given its order.
	 */
	inline std::size_t coupleAtOrd(std::size_t ord) const {
		NodeIdType i = std::floor(
				(2 * nbNodes - 1
						- std::sqrt(
								4 * nbNodes * nbNodes - 4 * nbNodes - 8 * ord
										+ 1)) / 2);
		std::size_t t = i * (2 * nbNodes - i - 3) / 2;
		NodeIdType j = ord - t + 1;
		return makeEdge(i, j);
	}

	/**
	 *@param iid An internal node ID.
	 *@return An iterator to one past the last neighbor of iid.
	 */
	EdgeIterator neighborsEnd(NodeIdType const &iid) const {
		return sja.begin() + sia[iid + 1];
	}

	/**
	 * @param iid The node internal ID.
	 * @return The degree of node iid.
	 */
	std::size_t getDeg(NodeIdType const &iid) const {
		return sia[iid + 1] - sia[iid];
	}

	/**
	 * Check if an edge exists in O(k_max).
	 * @param ii An internal node ID.
	 * @param ij An internal node ID.
	 * @return true if the edge (ii, ij) exists, false otherwise.
	 */
	bool isEdge(NodeIdType const &ii, NodeIdType const &ij) const {
		return isEdge(makeEdge(ii, ij));
	}

	/**
	 * Check if an edge exists in O(k_max).
	 * @param edge An edge.
	 * @return True if edge exists in the network, false otherwise.
	 */
	bool isEdge(EdgeType const &edge) const;

	/**
	 * @return The number of nodes in the network.
	 */
	std::size_t getNbNodes() const {
		return nbNodes;
	}

	/**
	 * @return The number of couples in the network.
	 */
	std::size_t getNbCouples() const {
		return nbNodes * (nbNodes - 1);
	}

	/**
	 * @return The number of edges in the network.
	 */
	std::size_t getNbEdges() const {
		return nbEdges;
	}

	/**
	 * @return The number of non-edges in the network.
	 */
	std::size_t getNbNonEdges() const {
		return nbNonEdges;
	}

	/**
	 * @return Average degree. Can only be called after the network is assembled.
	 */
	double getAvgDeg() const {
		return avgDeg;
	}

	/**
	 * @return Maximum degree. Can only be called after the network is assembled.
	 */
	std::size_t getMaxDeg() const {
		return maxDeg;
	}

	/**
	 * @return Minimum degree. Can only be called after the network is assembled.
	 */
	std::size_t getMinDeg() const {
		return minDeg;
	}

	/**
	 * Assemble the network. No changes to the network are allowed after calling this method.
	 */
	void assemble();

	/**
	 * Shuffles the nodes' internal IDs. This is useful to eliminate bias in methods that depend
	 * on node/edge order. Upon calling this method, all iterators and maps associated with the
	 * network are invalidated.
	 * @param seed Random number generator's seed.
	 */
	void shuffle(long int seed);

	/**
	 * @return A read-only (constant) external node ID iterator that points to the first node.
	 */
	LabelIterator labelsBegin() const {
		return nodesIdMap.kbegin();
	}

	/**
	 * @return A read-only (constant) external node ID iterator that points one past the last node.
	 */
	LabelIterator labelsEnd() const {
		return nodesIdMap.kend();
	}

	/**
	 * @return A read-only (constant) external node ID iterator that points to the first node.
	 */
	NodeIterator nodesBegin() const {
		return nodesIdMap.pbegin();
	}

	/**
	 * @return A read-only (constant) external node ID iterator that points one past the last node.
	 */
	NodeIterator nodesEnd() const {
		return nodesIdMap.pend();
	}

	/**
	 * @return A read-only (constant) external node ID iterator that points to the first node-degree couple.
	 */
	NodeDegIterator nodesDegBegin() const {
		return NodeDegIterator(*this, nodesBegin());
	}

	/**
	 * @return A read-only (constant) external node ID iterator that points one past the last node-degree couple.
	 */
	NodeDegIterator nodesDegEnd() const {
		return NodeDegIterator(*this, nodesEnd());
	}

	/**
	 * @return A read-only (constant) iterator that points to the first edge (with internal ID).
	 */
	EdgeIterator edgesBegin() const {
		return aja.begin();
	}

	/**
	 * @return A read-only (constant) iterator that points one past the last edge (with internal ID).
	 */
	EdgeIterator edgesEnd() const {
		return aja.end();
	}

	/**
	 * @return A read-only (constant) iterator that points to the first non-edge (with internal ID).
	 */
	NonEdgeIterator nonEdgesBegin() const {
		return NonEdgeIterator(*this, 0);
	}

	/**
	 * @return A read-only (constant) iterator that points one past the last non-edge (with internal ID).
	 */
	NonEdgeIterator nonEdgesEnd() const {
		return NonEdgeIterator(*this, nbNonEdges);
	}

	/**
	 * @param ratio Ratio of nodes that are selected.
	 * @param seed The random number gnerator's seed.
	 * @return a read-only (constant) randomized iterator that points to the first node (with internal ID).
	 */
	RndNodeIterator rndNodesBegin(double ratio, long int seed) const {
		return RndNodeIterator(nodesIdMap.pbegin(), nodesIdMap.pend(), ratio,
				seed);
	}

	/**
	 * @return A read-only (constant) randomized iterator that points one past the last node (with internal ID).
	 */
	RndNodeIterator rndNodesEnd() const {
		return RndNodeIterator(nodesIdMap.pend(), nodesIdMap.pend(), 0, 0);
	}

	/**
	 * @param ratio Ratio of nonedges that are selected.
	 * @param seed The random number gnerator's seed.
	 * @return a read-only (constant) randomized iterator that points to the first non-edge (with internal ID).
	 */
	RndNonEdgeIterator rndNonEdgesBegin(double ratio, long int seed) const {
		return RndNonEdgeIterator(*this, 0, ratio, seed);
	}

	/**
	 * @return A read-only (constant) randomized iterator that points one past the last non-edge (with internal ID).
	 */
	RndNonEdgeIterator rndNonEdgesEnd() const {
		return RndNonEdgeIterator(*this, nbNonEdges);
	}

	/**
	 * @param ratio Ratio of edges that are selected.
	 * @param seed The random number gnerator's seed.
	 * @return a read-only (constant) randomized iterator that points to the first edge (with internal ID).
	 */
	RndEdgeIterator rndEdgesBegin(double ratio, long int seed) const {
		return RndEdgeIterator(aja.begin(), aja.end(), ratio, seed);
	}

	/**
	 * @return A read-only (constant) randomized iterator that points one past the last edge (with internal ID).
	 */
	RndEdgeIterator rndEdgesEnd() const {
		return RndEdgeIterator(aja.end(), aja.end(), 0, 0);
	}

	/**
	 * Compute some degree statistics.
	 * @param minDeg (output parameter) minimum degree.
	 * @param maxDeg (output parameter) maximum degree.
	 * @param avgDeg (output parameter) average degree.
	 */
	void getDegStat(std::size_t &minDeg, std::size_t &maxDeg,
			double &avgDeg) const {
		minDeg = this->minDeg;
		maxDeg = this->maxDeg;
		avgDeg = this->avgDeg;
	}

	/**
	 * @tparam ValueT Value type.
	 * @return A node map.
	 */
	template<typename ValueT> NodeMap<ValueT> createNodeMap() const {
		return NodeMap<ValueT>(*this);
	}

	/**
	 * @tparam ValueT Value type.
	 * @return A pointer to a node map.
	 */
	template<typename ValueT> NodeMapSP<ValueT> createNodeMapSP() const {
		return std::shared_ptr<NodeMap<ValueT>>(new NodeMap<ValueT>(*this));
	}

	/**
	 * @tparam ValueT Value type.
	 * @return A sparse node map.
	 */
	template<typename ValueT> NodeSMap<ValueT> createNodeSMap(
			ValueT const &defVal) const {
		return NodeSMap<ValueT>(*this, defVal);
	}

	/**
	 * @tparam ValueT Value type.
	 * @return A pointer to a sparse node map.
	 */
	template<typename ValueT> NodeSMapSP<ValueT> createNodeSMapSP(
			ValueT const &defVal) const {
		return std::shared_ptr<NodeSMap<ValueT>>(
				new NodeSMap<ValueT>(*this, defVal));
	}

	/**
	 * @tparam ValueT Value type.
	 * @return An edge map.
	 */
	template<typename ValueT> EdgeMap<ValueT> createEdgeMap() const {
		return EdgeMap<ValueT>(*this);
	}

	/**
	 * @tparam ValueT Value type.
	 * @return A pointer to an edge map.
	 */
	template<typename ValueT> EdgeMapSP<ValueT> createEdgeMapSP() const {
		return std::shared_ptr<EdgeMap<ValueT>>(new EdgeMap<ValueT>(*this));
	}

	/**
	 * @param i A node ID.
	 * @param j A node ID.
	 * @return The number of common neighbors between i and j.
	 */
	std::size_t getNbCommonNeighbors(NodeIdType const &i,
			NodeIdType const &j) const {
		std::vector<EdgeType> c;
		getCommonNeighbors(i, j, std::back_inserter(c));
		return c.size();
	}

	/**
	 * @param i A node ID.
	 * @param j A node ID.
	 * @param inserter An inserter iterator to insert the common neighbors of i and j.
	 */
	template<typename InserterIterator> void getCommonNeighbors(
			NodeIdType const &i, NodeIdType const &j,
			InserterIterator inserter) const {
		std::set_intersection(sja.begin() + sia[i], sja.begin() + sia[i + 1],
				sja.begin() + sia[j], sja.begin() + sia[j + 1], inserter,
				compareEdgeEnd);
	}

	/**
	 * @param i A node ID.
	 * @param j A node ID.
	 * @return The set of common neighbors of i and j as an std::set.
	 */
	std::set<NodeIdType> getCommonNeighbors(NodeIdType const &i,
			NodeIdType const &j) const {

		std::vector<EdgeType> c;
		std::set_intersection(sja.begin() + sia[i], sja.begin() + sia[i + 1],
				sja.begin() + sia[j], sja.begin() + sia[j + 1],
				std::back_inserter(c), compareEdgeEnd);
		std::set<NodeIdType> ns;
		for (auto it = c.begin(); it != c.end(); ++it) {
			ns.insert(end(*it));
		}
		return ns;
	}

	/**
	 * @param i A node ID.
	 * @return The clustering coefficient of i. 
	 */
	double getCC(NodeIdType const &i) const {

		std::size_t nbn = getDeg(i);
		if (nbn <= 1) {
			return 0;
		}
		std::size_t nbc = 0;
		for (auto fit = neighborsBegin(i); fit < neighborsEnd(i); ++fit) {
			auto j = end(*fit);
			for (auto sit = (fit + 1); sit < neighborsEnd(i); ++sit) {
				auto k = end(*sit);
				if (isEdge(j, k)) {
					nbc++;
				}
			}
		}
		return (2.0 * nbc) / (nbn * (nbn - 1.0));
	}

	/**
	 * @return The average clustering coefficient of the network. 
	 */
	double getCC() const {
		double cc = 0;
		for (auto it = nodesBegin(); it != nodesEnd(); ++it) {
			cc += getCC(it->first);
		}
		return cc / nbNodes;
	}

	/**
	 * @param srcId The source node ID.
	 * @param endId The end node ID.
	 * @param length Specified length.
	 * @return The number of paths of length exactly length joining srcNode and endNode.
	 */
	std::size_t getNbPaths(NodeIdType const &srcId, NodeIdType const &endId,
			std::size_t length) const;

	/**
	 * @param ns A set of nodes.
	 * @return The number of edges connecting nodes in the set ns.
	 */
	std::size_t getNbInEdges(std::set<NodeIdType> const &ns) const;

	/**
	 * Read network from file.
	 * @param fileName The file name.
	 * @param ignoreRepetitions Whether to ignore repeated edges.
	 * @param ignoreLoops Whether to ignore loops.
	 * @return The read network.
	 */
	static std::shared_ptr<UNetwork<LabelType, NodeIdType, EdgeType, Comparator>> read(
			std::string fileName, bool ignoreRepetitions = false,
			bool ignoreLoops = false);

	/**
	 * Read couples from file.
	 * @param fileName The file name.
	 * @return The edges.
	 */
	std::shared_ptr<std::vector<EdgeType>> readCouples(
			std::string fileName) const;

	/**
	 * Generate an Erdos-Renyi (random) network.
	 * 	@param nbNodes The number of nodes.
	 * 	@param pr The probability of connecting any two nodes.
	 * 	@param seed Number generator seed.
	 */
	static std::shared_ptr<
			UNetwork<unsigned int, NodeIdType, EdgeType, std::less<unsigned int>>> generateERN(
			std::size_t nbNodes, double pr, long int seed);

	/**
	 * Generate an random connected network.
	 * 	@param nbNodes The number of nodes.
	 * 	@param pr The probability of connecting any two nodes.
	 * 	@param seed Number generator seed.
	 */
	static std::shared_ptr<
			UNetwork<unsigned int, NodeIdType, EdgeType, std::less<unsigned int>>> generateRNC(
			std::size_t nbNodes, double pr, long int seed);

	/**
	 * Generate a regular two-dimensional grid of size \f$ h \times w \f$, where \f$ h\f$ is the height of grid
	 * and \f$ w\f$ its width. The nodes are connected to their four neighbors right, left, up and down, except
	 * of course for boundary nodes. The nodes are assigned coordinates in \f$ [0,1]^2\f$ using equi-distant
	 * spacing in each dimension. The step size in the first dimension is \f$1/(h-1)\f$, in the second dimension,
	 * it is \f$1/(w-1)\f$.
	 * 	@param h Height of the grid.
	 * 	@param w Width of the grid.
	 */
	static std::shared_ptr<
			UNetwork<unsigned int, NodeIdType, EdgeType, std::less<unsigned int>>> generateREG(
			std::size_t h, std::size_t w);

	/**
	 * Write adjacency matrix in sparse form to file.
	 * @param fileName the file name.
	 */
	void write(std::string fileName) const;

	/**
	 * Print edges to std::cout.
	 */
	void print() const;

	/**
	 * Print edges to std::cout.
	 * @param edgesBegin Iterator to the beginning of the edges.
	 * @param edgesEnd Iterator to one past the the end of the edges.
	 */
	template<typename ForwardIterator> void printEdges(
			ForwardIterator edgesBegin, ForwardIterator edgesEnd) const {
		for (auto it = edgesBegin; it != edgesEnd; ++it) {
			std::cout << getLabel(start(*it)) << "\t" << getLabel(end(*it))
					<< std::endl;
		}
	}

	/**
	 * Destructor.
	 */
	virtual ~UNetwork() = default;
}
;

}
/* namespace LinkPred */

#endif /* UNETWORK_HPP_ */