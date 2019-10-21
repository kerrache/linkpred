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
 * @brief Contains the implementation of test data related classes.
 */

#ifndef NETWORKMANIPULATOR_HPP_
#define NETWORKMANIPULATOR_HPP_

#include "LinkPredConfig.hpp"
#include "linkpred/core/unetwork.hpp"
#include "linkpred/core/dnetwork.hpp"
#include "linkpred/utils/utilities.hpp"
#include "linkpred/core/graphtraversal.hpp"
#include "linkpred/utils/log.hpp"
#include <memory>
#include <vector>
#include <set>
#include <iterator>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

namespace LinkPred {

/**
 * @brief Enumeration of all classes of links.
 */
enum LinkClass {
	TP, /**< True positive link. */
	FN, /**< False negative link. */
	FP, /**< False positive link. */
	TN /**< True negative link. */
};

/**
 * @brief Generate true positives and true negatives.
 * @tparam NetworkT Network type.
 * @tparam EdgeContainerT The container used to store edges.
 */
template<typename NetworkT = UNetwork<>, typename EdgeContainerT = std::vector<
		typename NetworkT::EdgeType>> class TestEdgeGenerator {
	using NetworkSP = std::shared_ptr<NetworkT>; /**< Shared pointer to a network. */
	using NetworkCSP = std::shared_ptr<const NetworkT>; /**< Constant shared pointer to a network. */
	using EdgeType = typename NetworkT::EdgeType; /**< The edge type. */
	using NonEdgeIterator = typename NetworkT::NonEdgeIterator; /**< Negative edge iterator type. */
	using EdgeIterator = typename NetworkT::EdgeIterator; /**< Edge iterator type. */

protected:
	NetworkCSP refNet; /**< Reference network. */
	NetworkCSP obsNet; /**< Observed network. */
	std::set<EdgeType> addLinksSet; /**< Set of added links. */
	bool aTP = false; /**< Use all TP. */
	double tpRatio = 0; /**< TP ratio. */
	bool aTN = false; /**< Use all TN. */
	double tnRatio = 0; /**< TN ratio. */
	long int tpSeed; /**< Seed for TP. */
	long int tnSeed; /**< Seed for TN. */
	typename NetworkT::NonEdgeIterator tnEndIt; /**< End of TN edges. Cached for performance purposes. */
	typename NetworkT::RndNonEdgeIterator tnRndEndIt; /**< End of TN random edges. Cached for performance purposes. */
	typename NetworkT::EdgeIterator tpEndIt; /**< End of TP edges. Cached for performance purposes. */
	typename NetworkT::RndEdgeIterator tpRndEndIt; /**< End of TP random edges. Cached for performance purposes. */

public:

	/**
	 * @brief TP edges iterator.
	 */
	class TPEdgeIterator: public std::iterator<std::random_access_iterator_tag,
			const EdgeType, long int> {

		friend class TestEdgeGenerator; /**< TestEdgeGenerator is a friend. */

	protected:
		TestEdgeGenerator const &tg; /**< Test generator object (owner). */
		typename NetworkT::EdgeIterator eit; /**< Edge iterator. */
		typename NetworkT::RndEdgeIterator reit; /**< Random edge iterator. */

		/**
		 * Constructor.
		 * @param tg The test generator object (owner).
		 */
		TPEdgeIterator(TestEdgeGenerator const &tg) :
				tg(tg), eit(tg.obsNet->edgesBegin()), reit(
						tg.obsNet->rndEdgesBegin(tg.tpRatio, tg.tpSeed)) {
			if (tg.aTP) {
				while ((eit != tg.obsNet->edgesEnd())
						&& ((tg.addLinksSet.count(*eit) != 0)
								|| (tg.addLinksSet.count(
										NetworkT::reverseEdge(*eit)) != 0))) {
					++eit;
				}
			} else {
				while ((reit != tg.obsNet->rndEdgesEnd())
						&& ((tg.addLinksSet.count(*reit) != 0)
								|| (tg.addLinksSet.count(
										NetworkT::reverseEdge(*reit)) != 0))) {
					++reit;
				}
			}
		}

		/**
		 * Constructor.
		 * @param tg The test generator object (owner).
		 */
		TPEdgeIterator(TestEdgeGenerator const &tg,
				typename NetworkT::EdgeIterator const &eit,
				typename NetworkT::RndEdgeIterator const &reit) :
				tg(tg), eit(eit), reit(reit) {
		}

	public:
		using pointer = typename std::iterator<std::random_access_iterator_tag, const EdgeType, long int>::pointer; /**< The pointer type associated with the iterator. */
		using reference = typename std::iterator<std::random_access_iterator_tag, const EdgeType, long int>::reference; /**< The reference type associated with the iterator. */
		using difference_type = typename std::iterator<std::random_access_iterator_tag, const EdgeType, long int>::difference_type; /**< The difference type associated with the iterator. */

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		TPEdgeIterator(TPEdgeIterator const &that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		TPEdgeIterator& operator =(TPEdgeIterator const &that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		TPEdgeIterator(TPEdgeIterator &&that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		TPEdgeIterator& operator =(TPEdgeIterator &&that) = default;

		/**
		 * Dereference operator.
		 * @return A reference to the object to which the iterator points.
		 */
		reference operator*() {
			if (tg.aTP) {
				return *eit;
			} else {
				return *reit;
			}
		}

		/**
		 * Arrow operator.
		 * @return A pointer to the object to which the iterator points.
		 */
		pointer operator->() {
			if (tg.aTP) {
				return &(*eit);
			} else {
				return &(*reit);
			}
		}

		/**
		 * Pre-increment operator.
		 * @return A reference to the new iterator.
		 */
		TPEdgeIterator& operator++() {
			if (tg.aTP) {
				++eit;
				while ((eit != tg.obsNet->edgesEnd())
						&& ((tg.addLinksSet.count(*eit) != 0)
								|| (tg.addLinksSet.count(
										NetworkT::reverseEdge(*eit)) != 0))) {
					++eit;
				}
			} else {
				++reit;
				while ((reit != tg.obsNet->rndEdgesEnd())
						&& ((tg.addLinksSet.count(*reit) != 0)
								|| (tg.addLinksSet.count(
										NetworkT::reverseEdge(*reit)) != 0))) {
					++reit;
				}
			}
			return *this;
		}

		/**
		 * Pre-decrement operator.
		 * @return A reference to the new iterator.
		 */
		TPEdgeIterator& operator--() {
			if (tg.aTP) {
				if (tg.addLinksSet.size() > 0) {
					throw std::runtime_error(
							"Operator -- not supported when the set of added links is not empty.");
				}
				--eit;
			} else {
				throw std::runtime_error(
						"Operator -- not supported for randomized iterators.");
			}
			return *this;
		}

		/**
		 * Post-increment operator.
		 * @return A reference to the new iterator.
		 */
		TPEdgeIterator operator++(int) {
			auto that = *this;
			++(*this);
			return that;
		}

		/**
		 * Post-decrement operator.
		 * @return A reference to the new iterator.
		 */
		TPEdgeIterator operator--(int) {
			auto that = *this;
			--(*this);
			return that;
		}

		/**
		 * Arithmetic + operator.
		 * @param n Increment value.
		 * @return The new iterator.
		 */
		TPEdgeIterator operator+(const difference_type &n) const {
			if (tg.aTP) {
				if (tg.addLinksSet.size() > 0) {
					throw std::runtime_error(
							"Operator +n not supported when the set of added links is not empty.");
				}
				auto that = *this;
				that.eit += n;
				return that;
			} else {
				throw std::runtime_error(
						"Operator +n not supported for randomized iterators.");
			}
		}

		/**
		 * Arithmetic += operator.
		 * @param n Increment value.
		 * @return A reference to the new iterator.
		 */
		TPEdgeIterator& operator+=(const difference_type &n) {
			if (tg.aTP) {
				if (tg.addLinksSet.size() > 0) {
					throw std::runtime_error(
							"Operator +=n not supported when the set of added links is not empty.");
				}
				eit += n;
				return *this;
			} else {
				throw std::runtime_error(
						"Operator +=n not supported for randomized iterators.");
			}
		}

		/**
		 * Arithmetic - operator.
		 * @param n Decrement value.
		 * @return The new iterator.
		 */
		TPEdgeIterator operator-(const difference_type &n) const {
			if (tg.aTP) {
				if (tg.addLinksSet.size() > 0) {
					throw std::runtime_error(
							"Operator -n not supported when the set of added links is not empty.");
				}
				auto that = *this;
				that.eit -= n;
				return that;
			} else {
				throw std::runtime_error(
						"Operator -n not supported for randomized iterators.");
			}
		}

		/**
		 * Arithmetic -= operator.
		 * @param n Decrement value.
		 * @return A reference to the new iterator.
		 */
		TPEdgeIterator& operator-=(const difference_type &n) {
			if (tg.aTP) {
				if (tg.addLinksSet.size() > 0) {
					throw std::runtime_error(
							"Operator -=n not supported when the set of added links is not empty.");
				}
				eit -= n;
				return *this;
			} else {
				throw std::runtime_error(
						"Operator -=n not supported for randomized iterators.");
			}
		}

		/**
		 * Difference between the present iterator and the one passed as parameter.
		 * @param that The other iterator.
		 * @return The difference between the current and that iterator.
		 */
		difference_type operator-(const TPEdgeIterator &that) const {
			if (tg.aTP) {
				if (tg.addLinksSet.size() > 0) {
					throw std::runtime_error(
							"Operator - not supported when the set of added links is not empty.");
				}
				return eit - that.eit;
			} else {
				throw std::runtime_error(
						"Operator - not supported for randomized iterators.");
			}
		}

		/**
		 * @param that The other iterator.
		 * @return True if this equals that.
		 */
		bool operator==(const TPEdgeIterator &that) const {
			if (tg.aTP != that.tg.aTP) {
				return false;
			}
			if (tg.aTP) {
				return that.eit == eit;
			} else {
				return that.reit == reit;
			}
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is not equal to that.
		 */
		bool operator!=(const TPEdgeIterator &that) const {
			return !(*this == that);
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is less than that.
		 */
		bool operator<(const TPEdgeIterator &that) const {
			if (tg.aTP != that.tg.aTP) {
				return false;
			}
			if (tg.aTP) {
				return eit < that.eit;
			} else {
				throw std::runtime_error(
						"Operator < not supported for randomized iterators.");
			}
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is greater that.
		 */
		bool operator>(const TPEdgeIterator &that) const {
			if (tg.aTP != that.tg.aTP) {
				return false;
			}
			if (tg.aTP) {
				return eit > that.eit;
			} else {
				throw std::runtime_error(
						"Operator > not supported for randomized iterators.");
			}
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is greater or equal to that.
		 */
		bool operator<=(const TPEdgeIterator &that) const {
			return !(*this > that);
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is less or equal to that.
		 */
		bool operator>=(const TPEdgeIterator &that) const {
			return !(*this < that);
		}
	};

	/**
	 * @brief TN edges iterator.
	 */
	class TNEdgeIterator: public std::iterator<std::random_access_iterator_tag,
			const EdgeType, long int> {

		friend class TestEdgeGenerator; /**< TestEdgeGenerator is a friend. */

	protected:
		TestEdgeGenerator const &tg; /**< Test generator object (owner). */
		bool aTN = true; /**< Use all true negatives? */
		NetworkCSP refNet; /**< Reference network. */
		typename NetworkT::NonEdgeIterator neit; /**< Non-edge iterator. */
		typename NetworkT::RndNonEdgeIterator rneit; /**< Random non-edge iterator. */

		/**
		 * Constructor.
		 * @param tg The test generator object (owner).
		 */
		TNEdgeIterator(TestEdgeGenerator const &tg,
				typename NetworkT::NonEdgeIterator const &_neit,
				typename NetworkT::RndNonEdgeIterator const &_rneit) :
				tg(tg), aTN(tg.aTN), refNet(tg.refNet), neit(_neit), rneit(
						_rneit) {
			if (aTN) {
				while ((neit != tg.refNet->nonEdgesEnd())
						&& ((tg.addLinksSet.count(*neit) != 0)
								|| (tg.addLinksSet.count(
										NetworkT::reverseEdge(*neit)) != 0))) {
					++neit;
				}
			} else {
				while ((rneit != tg.refNet->rndNonEdgesEnd())
						&& ((tg.addLinksSet.count(*rneit) != 0)
								|| (tg.addLinksSet.count(
										NetworkT::reverseEdge(*rneit)) != 0))) {
					++rneit;
				}
			}
		}

	public:
		using pointer = typename std::iterator<std::random_access_iterator_tag, const EdgeType, long int>::pointer; /**< The pointer type associated with the iterator. */
		using reference = typename std::iterator<std::random_access_iterator_tag, const EdgeType, long int>::reference; /**< The reference type associated with the iterator. */
		using difference_type = typename std::iterator<std::random_access_iterator_tag, const EdgeType, long int>::difference_type; /**< The difference type associated with the iterator. */

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		TNEdgeIterator(TNEdgeIterator const &that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		TNEdgeIterator& operator =(TNEdgeIterator const &that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		TNEdgeIterator(TNEdgeIterator &&that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		TNEdgeIterator& operator =(TNEdgeIterator &&that) = default;

		/**
		 * Dereference operator.
		 * @return A reference to the object to which the iterator points.
		 */
		reference operator*() {
			if (aTN) {
				return *neit;
			} else {
				return *rneit;
			}
		}

		/**
		 * Arrow operator.
		 * @return A pointer to the object to which the iterator points.
		 */
		pointer operator->() {
			if (aTN) {
				return &(*neit);
			} else {
				return &(*rneit);
			}
		}

		/**
		 * Pre-increment operator.
		 * @return A reference to the new iterator.
		 */
		TNEdgeIterator& operator++() {
			EdgeType e;
			if (aTN) {
				do {
					++neit;
					e = *neit;
				} while ((tg.addLinksSet.count(e) != 0)
						|| (tg.addLinksSet.count(NetworkT::reverseEdge(e)) != 0));
			} else {
				do {
					++rneit;
					e = *rneit;
				} while ((tg.addLinksSet.count(e) != 0)
						|| (tg.addLinksSet.count(NetworkT::reverseEdge(e)) != 0));
			}
			return *this;
		}

		/**
		 * Pre-decrement operator.
		 * @return A reference to the new iterator.
		 */
		TNEdgeIterator& operator--() {
			if (aTN) {
				if (tg.addLinksSet.size() > 0) {
					throw std::runtime_error(
							"Operator -- not supported when the set of added links is not empty.");
				}
				--neit;
			} else {
				throw std::runtime_error(
						"Operator -- not supported for randomized iterators.");
			}
			return *this;
		}

		/**
		 * Post-increment operator.
		 * @return A reference to the new iterator.
		 */
		TNEdgeIterator operator++(int) {
			auto that = *this;
			++(*this);
			return that;
		}

		/**
		 * Post-decrement operator.
		 * @return A reference to the new iterator.
		 */
		TNEdgeIterator operator--(int) {
			auto that = *this;
			--(*this);
			return that;
		}

		/**
		 * Arithmetic + operator.
		 * @param n Increment value.
		 * @return The new iterator.
		 */
		TNEdgeIterator operator+(const difference_type &n) const {
			if (aTN) {
				if (tg.addLinksSet.size() > 0) {
					throw std::runtime_error(
							"Operator +n not supported when the set of added links is not empty.");
				}
				auto that = *this;
				that.neit += n;
				return that;
			} else {
				throw std::runtime_error(
						"Operator +n not supported for randomized iterators.");
			}
		}

		/**
		 * Arithmetic += operator.
		 * @param n Increment value.
		 * @return A reference to the new iterator.
		 */
		TNEdgeIterator& operator+=(const difference_type &n) {
			if (aTN) {
				if (tg.addLinksSet.size() > 0) {
					throw std::runtime_error(
							"Operator +=n not supported when the set of added links is not empty.");
				}
				neit += n;
				return *this;
			} else {
				throw std::runtime_error(
						"Operator +=n not supported for randomized iterators.");
			}
		}

		/**
		 * Arithmetic - operator.
		 * @param n Decrement value.
		 * @return The new iterator.
		 */
		TNEdgeIterator operator-(const difference_type &n) const {
			if (aTN) {
				if (tg.addLinksSet.size() > 0) {
					throw std::runtime_error(
							"Operator -n not supported when the set of added links is not empty.");
				}
				auto that = *this;
				that.neit -= n;
				return that;
			} else {
				throw std::runtime_error(
						"Operator -n not supported for randomized iterators.");
			}
		}

		/**
		 * Arithmetic -= operator.
		 * @param n Decrement value.
		 * @return A reference to the new iterator.
		 */
		TNEdgeIterator& operator-=(const difference_type &n) {
			if (aTN) {
				if (tg.addLinksSet.size() > 0) {
					throw std::runtime_error(
							"Operator -=n not supported when the set of added links is not empty.");
				}
				neit -= n;
				return *this;
			} else {
				throw std::runtime_error(
						"Operator -=n not supported for randomized iterators.");
			}
		}

		/**
		 * Difference between the present iterator and the one passed as parameter.
		 * @param that The other iterator.
		 * @return The difference between the current and that iterator.
		 */
		difference_type operator-(const TNEdgeIterator &that) const {
			if (aTN) {
				if (tg.addLinksSet.size() > 0) {
					throw std::runtime_error(
							"Operator - not supported when the set of added links is not empty.");
				}
				return neit - that.neit;
			} else {
				throw std::runtime_error(
						"Operator - not supported for randomized iterators.");
			}
		}

		/**
		 * @param that The other iterator.
		 * @return True if this equals that.
		 */
		bool operator==(const TNEdgeIterator &that) const {
			if (aTN != that.aTN) {
				return false;
			}
			if (aTN) {
				return that.neit == neit;
			} else {
				return that.rneit == rneit;
			}
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is not equal to that.
		 */
		bool operator!=(const TNEdgeIterator &that) const {
			return !(*this == that);
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is less than that.
		 */
		bool operator<(const TNEdgeIterator &that) const {
			if (aTN != that.aTN) {
				return false;
			}
			if (aTN) {
				return neit < that.neit;
			} else {
				throw std::runtime_error(
						"Operator < not supported for randomized iterators.");
			}
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is greater that.
		 */
		bool operator>(const TNEdgeIterator &that) const {
			if (aTN != that.aTN) {
				return false;
			}
			if (aTN) {
				return neit > that.neit;
			} else {
				throw std::runtime_error(
						"Operator > not supported for randomized iterators.");
			}
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is greater or equal to that.
		 */
		bool operator<=(const TNEdgeIterator &that) const {
			return !(*this > that);
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is less or equal to that.
		 */
		bool operator>=(const TNEdgeIterator &that) const {
			return !(*this < that);
		}
	};

	/**
	 * @param refNet The reference network.
	 * @param obsNet The network.
	 * @param remLinks Removed edges.
	 * @param addLinks Added edges.
	 * @param aTP Whether to use all true positive links in the test set.
	 * @param tpRatio Ratio of true positive inks to e used in the test set. This parameter is only relevant when aTP is false.
	 * @param aTN Whether to use all true negative links in the test set.
	 * @param tnRatio Ratio of true negative inks to e used in the test set. This parameter is only relevant when aTN is false.
	 * @param seed The random number generator's seed.
	 */
	TestEdgeGenerator(NetworkCSP refNet, NetworkCSP obsNet,
			std::shared_ptr<EdgeContainerT> remLinks,
			std::shared_ptr<EdgeContainerT> addLinks, bool aTP, double tpRatio,
			bool aTN, double tnRatio, long int seed) :
			refNet(refNet), obsNet(obsNet), aTP(aTP), tpRatio(tpRatio), aTN(
					aTN), tnRatio(tnRatio), tnEndIt(refNet->nonEdgesEnd()), tnRndEndIt(
					refNet->rndNonEdgesEnd()), tpEndIt(obsNet->edgesEnd()), tpRndEndIt(
					obsNet->rndEdgesEnd()) {
		if ((tpRatio < 0) || (tpRatio > 1)) {
			logger(logError, "tpRatio: " << tpRatio)
			throw std::domain_error("The ratios must be between 0 and 1");
		}

		if ((tnRatio < 0) || (tnRatio > 1)) {
			logger(logError, "tnRatio: " << tnRatio)
			throw std::domain_error("The ratios must be between 0 and 1");
		}

		for (auto it = addLinks->begin(); it != addLinks->end(); ++it) {
			addLinksSet.insert(*it);
		}
		RandomGen rng(seed);
		tpSeed = rng.getInt();
		tnSeed = rng.getInt();
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	TestEdgeGenerator(TestEdgeGenerator const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	TestEdgeGenerator& operator =(TestEdgeGenerator const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	TestEdgeGenerator(TestEdgeGenerator &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	TestEdgeGenerator& operator =(TestEdgeGenerator &&that) = default;

	/**
	 * Generate true positive links.
	 * @tparam OutputEdgeIteratorT Iterator write edges.
	 * @param oit Output edge iterator.
	 */
	template<typename OutputEdgeIteratorT = typename std::vector<
			typename NetworkT::EdgeType>::iterator> void generateTP(
			OutputEdgeIteratorT oit) {

		// Add + Remove
		// The true positive links are the edges of obsNet that were not added
		if (aTP) {
			for (auto it = obsNet->edgesBegin(); it != obsNet->edgesEnd();
					++it) {
				if ((addLinksSet.count(*it) == 0)
						&& (addLinksSet.count(NetworkT::reverseEdge(*it)) == 0)) {
					oit = *it;
					++oit;
				}
			}
		} else {
			for (auto it = obsNet->rndEdgesBegin(tpRatio, tpSeed);
					it != obsNet->rndEdgesEnd(); ++it) {
				if ((addLinksSet.count(*it) == 0)
						&& (addLinksSet.count(NetworkT::reverseEdge(*it)) == 0)) {
					oit = *it;
					++oit;
				}
			}
		}
	}

	/**
	 * Generate true negative links.
	 * @tparam OutputEdgeIteratorT Iterator write edges.
	 * @param oit Output edge iterator.
	 */
	template<typename OutputEdgeIteratorT = typename std::vector<
			typename NetworkT::EdgeType>::iterator> void generateTN(
			OutputEdgeIteratorT oit) {

		// Add + Remove
		// The true negative links are the non-edges of refNet that were not added
		if (aTN) {
			for (auto it = refNet->nonEdgesBegin(); it != refNet->nonEdgesEnd();
					++it) {
				if ((addLinksSet.count(*it) == 0)
						&& (addLinksSet.count(NetworkT::reverseEdge(*it)) == 0)) {
					oit = *it;
					++oit;
				}
			}
		} else {
			for (auto it = refNet->rndNonEdgesBegin(tnRatio, tnSeed);
					it != refNet->rndNonEdgesEnd(); ++it) {
				if ((addLinksSet.count(*it) == 0)
						&& (addLinksSet.count(NetworkT::reverseEdge(*it)) == 0)) {
					oit = *it;
					++oit;
				}
			}
		}
	}

	/**
	 * @return Iterator to the first true positive link.
	 */
	auto tpBegin() const {
		return TPEdgeIterator(*this);
	}

	/**
	 * @return Iterator to one-past the last true positive link.
	 */
	auto tpEnd() const {
		return TPEdgeIterator(*this, tpEndIt, tpRndEndIt);
	}

	/**
	 * @return Iterator to the first true positive link.
	 */
	auto tnBegin() const {
		return TNEdgeIterator(*this, refNet->nonEdgesBegin(),
				refNet->rndNonEdgesBegin(tnRatio, tnSeed));
	}

	/**
	 * @return Iterator to one-past the last true positive link.
	 */
	auto tnEnd() const {
		return TNEdgeIterator(*this, tnEndIt, tnRndEndIt);
	}

	/**
	 * Destructor.
	 */
	virtual ~TestEdgeGenerator() = default;
};

/**
 * @brief Test data.
 * @tparam NetworkT The network type.
 * @tparam EdgeContainerT The container used to store edges.
 */
template<typename NetworkT = UNetwork<>, typename EdgeContainerT = std::vector<
		typename NetworkT::EdgeType>> class TestData {
public:
	using NetworkSP = std::shared_ptr<NetworkT>; /**< Shared pointer to a network. */
	using NetworkCSP = std::shared_ptr<const NetworkT>; /**< Constant shared pointer to anetwork. */
	using EdgeType = typename NetworkT::EdgeType; /**< The edge type. */
	using NonEdgeIterator = typename NetworkT::NonEdgeIterator; /**< Negative edge iterator type. */
	using EdgeIterator = typename NetworkT::EdgeIterator; /**< Edge iterator type. */

protected:
	NetworkCSP refNet; /**< Reference network. */
	NetworkCSP obsNet; /**< Observed network. */
	std::shared_ptr<EdgeContainerT> remLinks; /**< Removed edges. */
	std::shared_ptr<EdgeContainerT> addLinks; /**< Added edges. */
	std::shared_ptr<EdgeContainerT> tpLinks; /**< True positive links. */
	std::shared_ptr<EdgeContainerT> tnLinks; /**< True negative links. */
	std::shared_ptr<EdgeContainerT> pos; /**< Positive links. */
	std::shared_ptr<EdgeContainerT> neg; /**< Negative links. */
	std::shared_ptr<std::set<EdgeType>> remLinksMap; /**< Map containing removed links. */
	std::shared_ptr<TestEdgeGenerator<NetworkT, EdgeContainerT>> eg; /**< Edge generator. */
	LinkClass posClass; /**< The class of edges used as the positive instances in the test set. */
	LinkClass negClass; /**< The class of edges used as the negative instances in the test set. */
	bool tpGenerated = false; /**< True positives are generated and stored in memory. */
	bool tnGenerated = false; /**< True negatives are generated and stored in memory. */
	bool locked = false; /**< Whether the test data is locked. If so, no modification can be done on the object. */

public:

	/**
	 * @brief Enumeration of iterator types.
	 */
	enum IteratorType {
		ECEIT, /**< Edge container iterator. */
		TPEIT, /**< True negative edge iterator. */
		TNEIT /**< True positive edge iterator. */
	};

	/**
	 * @brief Test edges iterator.
	 */
	class TestEdgeIterator: public std::iterator<
			std::random_access_iterator_tag, const EdgeType, long int> {

		friend class TestData; /**< TestData is a friend. */

	protected:
		IteratorType itType; /**< Iterator type. */
		typename EdgeContainerT::const_iterator eceit; /**< Edge container iterator. */
		typename TestEdgeGenerator<NetworkT, EdgeContainerT>::TPEdgeIterator tpeit; /**< True positive edge iterator. */
		typename TestEdgeGenerator<NetworkT, EdgeContainerT>::TNEdgeIterator tneit; /**< True negative edge iterator. */

	public:
		using pointer = typename std::iterator<std::random_access_iterator_tag, const EdgeType, long int>::pointer; /**< The pointer type associated with the iterator. */
		using reference = typename std::iterator<std::random_access_iterator_tag, const EdgeType, long int>::reference; /**< The reference type associated with the iterator. */
		using difference_type = typename std::iterator<std::random_access_iterator_tag, const EdgeType, long int>::difference_type; /**< The difference type associated with the iterator. */

		/**
		 * Constructor.
		 * @param td The test data object (owner).
		 */
		TestEdgeIterator(IteratorType const &itType,
				typename EdgeContainerT::const_iterator const &eceit,
				typename TestEdgeGenerator<NetworkT, EdgeContainerT>::TPEdgeIterator const &tpeit,
				typename TestEdgeGenerator<NetworkT, EdgeContainerT>::TNEdgeIterator const &tneit) :
				itType(itType), eceit(eceit), tpeit(tpeit), tneit(tneit) {
		}

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		TestEdgeIterator(TestEdgeIterator const &that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		TestEdgeIterator& operator =(TestEdgeIterator const &that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		TestEdgeIterator(TestEdgeIterator &&that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		TestEdgeIterator& operator =(TestEdgeIterator &&that) = default;

		/**
		 * Dereference operator.
		 * @return A reference to the object to which the iterator points.
		 */
		reference operator*() {
			switch (itType) {
			case ECEIT:
				return *eceit;
			case TPEIT:
				return *tpeit;
			case TNEIT:
				return *tneit;
			default:
				throw std::runtime_error("Unknown iterator type");
			}
		}

		/**
		 * Arrow operator.
		 * @return A pointer to the object to which the iterator points.
		 */
		pointer operator->() {
			switch (itType) {
			case ECEIT:
				return &(*eceit);
			case TPEIT:
				return &(*tpeit);
			case TNEIT:
				return &(*tneit);
			default:
				throw std::runtime_error("Unknown iterator type");
			}
		}

		/**
		 * Pre-increment operator.
		 * @return A reference to the new iterator.
		 */
		TestEdgeIterator& operator++() {
			switch (itType) {
			case ECEIT:
				++eceit;
				break;
			case TPEIT:
				++tpeit;
				break;
			case TNEIT:
				++tneit;
				break;
			default:
				throw std::runtime_error("Unknown iterator type");
			}
			return *this;
		}

		/**
		 * Pre-decrement operator.
		 * @return A reference to the new iterator.
		 */
		TestEdgeIterator& operator--() {
			switch (itType) {
			case ECEIT:
				--eceit;
				break;
			case TPEIT:
				--tpeit;
				break;
			case TNEIT:
				--tneit;
				break;
			default:
				throw std::runtime_error("Unknown iterator type");
			}
			return *this;
		}

		/**
		 * Post-increment operator.
		 * @return A reference to the new iterator.
		 */
		TestEdgeIterator operator++(int) {
			auto that = *this;
			++(*this);
			return that;
		}

		/**
		 * Post-decrement operator.
		 * @return A reference to the new iterator.
		 */
		TestEdgeIterator operator--(int) {
			auto that = *this;
			--(*this);
			return that;
		}

		/**
		 * Arithmetic + operator.
		 * @param n Increment value.
		 * @return The new iterator.
		 */
		TestEdgeIterator operator+(const difference_type &n) const {
			auto that = *this;
			switch (itType) {
			case ECEIT:
				that.eceit += n;
				break;
			case TPEIT:
				that.tpeit += n;
				break;
			case TNEIT:
				that.tneit += n;
				break;
			default:
				throw std::runtime_error("Unknown iterator type");
			}
			return that;
		}

		/**
		 * Arithmetic += operator.
		 * @param n Increment value.
		 * @return A reference to the new iterator.
		 */
		TestEdgeIterator& operator+=(const difference_type &n) {
			switch (itType) {
			case ECEIT:
				eceit += n;
				break;
			case TPEIT:
				tpeit += n;
				break;
			case TNEIT:
				tneit += n;
				break;
			default:
				throw std::runtime_error("Unknown iterator type");
			}
			return *this;
		}

		/**
		 * Arithmetic - operator.
		 * @param n Decrement value.
		 * @return The new iterator.
		 */
		TestEdgeIterator operator-(const difference_type &n) const {
			auto that = *this;
			switch (itType) {
			case ECEIT:
				that.eceit -= n;
				break;
			case TPEIT:
				that.tpeit -= n;
				break;
			case TNEIT:
				that.tneit -= n;
				break;
			default:
				throw std::runtime_error("Unknown iterator type");
			}
			return that;
		}

		/**
		 * Arithmetic -= operator.
		 * @param n Decrement value.
		 * @return A reference to the new iterator.
		 */
		TestEdgeIterator& operator-=(const difference_type &n) {
			switch (itType) {
			case ECEIT:
				eceit -= n;
				break;
			case TPEIT:
				tpeit -= n;
				break;
			case TNEIT:
				tneit -= n;
				break;
			default:
				throw std::runtime_error("Unknown iterator type");
			}
			return *this;
		}

		/**
		 * Difference between the present iterator and the one passed as parameter.
		 * @param that The other iterator.
		 * @return The difference between the current and that iterator.
		 */
		difference_type operator-(const TestEdgeIterator &that) const {
			switch (itType) {
			case ECEIT:
				return eceit - that.eceit;
			case TPEIT:
				return tpeit - that.tpeit;
			case TNEIT:
				return tneit - that.tneit;
			default:
				throw std::runtime_error("Unknown iterator type");
			}
		}

		/**
		 * @param that The other iterator.
		 * @return True if this equals that.
		 */
		bool operator==(const TestEdgeIterator &that) const {
			if (itType != that.itType) {
				return false;
			}
			switch (itType) {
			case ECEIT:
				return eceit == that.eceit;
			case TPEIT:
				return tpeit == that.tpeit;
			case TNEIT:
				return tneit == that.tneit;
			default:
				throw std::runtime_error("Unknown iterator type");
			}
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is not equal to that.
		 */
		bool operator!=(const TestEdgeIterator &that) const {
			return !(*this == that);
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is less than that.
		 */
		bool operator<(const TestEdgeIterator &that) const {
			if (itType != that.itType) {
				return false;
			}
			switch (itType) {
			case ECEIT:
				return eceit < that.eceit;
			case TPEIT:
				return tpeit < that.tpeit;
			case TNEIT:
				return tneit < that.tneit;
			default:
				throw std::runtime_error("Unknown iterator type");
			}
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is greater that.
		 */
		bool operator>(const TestEdgeIterator &that) const {
			if (itType != that.itType) {
				return false;
			}
			switch (itType) {
			case ECEIT:
				return eceit > that.eceit;
			case TPEIT:
				return tpeit > that.tpeit;
			case TNEIT:
				return tneit > that.tneit;
			default:
				throw std::runtime_error("Unknown iterator type");
			}
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is greater or equal to that.
		 */
		bool operator<=(const TestEdgeIterator &that) const {
			return !(*this > that);
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is less or equal to that.
		 */
		bool operator>=(const TestEdgeIterator &that) const {
			return !(*this < that);
		}
	};

protected:
	std::shared_ptr<TestEdgeIterator> posStrmEndItSP; /**< Pointer to end of positive links. Cached for performance purposes. */
	std::shared_ptr<TestEdgeIterator> negStrmEndItSP; /**< Pointer to end of negative links. Cached for performance purposes. */

public:
	/**
	 * @param refNet The reference network.
	 * @param obsNet The network.
	 * @param remLinks Removed edges.
	 * @param addLinks Added edges.
	 * @param tpLinks True positive links.
	 * @param tnLinks True negative links.
	 * @param posClass The class of edges used as the positive instances in the test set.
	 * @param negClass The class of edges used as the negative instances in the test set.
	 */
	TestData(NetworkCSP refNet, NetworkCSP obsNet,
			std::shared_ptr<EdgeContainerT> remLinks,
			std::shared_ptr<EdgeContainerT> addLinks,
			std::shared_ptr<EdgeContainerT> tpLinks,
			std::shared_ptr<EdgeContainerT> tnLinks, LinkClass posClass,
			LinkClass negClass) :
			refNet(refNet), obsNet(obsNet), remLinks(remLinks), addLinks(
					addLinks), tpLinks(tpLinks), tnLinks(tnLinks), posClass(
					posClass), negClass(negClass) {

		remLinksMap = std::make_shared<std::set<EdgeType>>();
		remLinksMap->insert(remLinks->cbegin(), remLinks->cend());
		if (remLinksMap->size() != remLinks->size()) {
			throw std::runtime_error(
					"Could not insert all removed links (are there duplicates?)");
		}
		switch (posClass) {
		case TP:
			pos = tpLinks;
			break;
		case FN:
			pos = remLinks;
			break;
		case FP:
			pos = addLinks;
			break;
		case TN:
			pos = tnLinks;
			break;
		default:
			throw std::runtime_error("Unknown link class");
		}

		switch (negClass) {
		case TP:
			neg = tpLinks;
			break;
		case FN:
			neg = remLinks;
			break;
		case FP:
			neg = addLinks;
			break;
		case TN:
			neg = tnLinks;
			break;
		default:
			throw std::runtime_error("Unknown link class");
		}

		tpGenerated = true;
		tnGenerated = true;
	}

	/**
	 * @param refNet The reference network.
	 * @param obsNet The network.
	 * @param remLinks Removed edges.
	 * @param addLinks Added edges.
	 * @param eg Edge generator.
	 * @param posClass The class of edges used as the positive instances in the test set.
	 * @param negClass The class of edges used as the negative instances in the test set.
	 */
	TestData(NetworkCSP refNet, NetworkCSP obsNet,
			std::shared_ptr<EdgeContainerT> remLinks,
			std::shared_ptr<EdgeContainerT> addLinks,
			std::shared_ptr<TestEdgeGenerator<NetworkT, EdgeContainerT>> eg,
			LinkClass posClass, LinkClass negClass) :
			refNet(refNet), obsNet(obsNet), remLinks(remLinks), addLinks(
					addLinks), eg(eg), posClass(posClass), negClass(negClass) {

		remLinksMap = std::make_shared<std::set<EdgeType>>();
		remLinksMap->insert(remLinks->cbegin(), remLinks->cend());
		if (remLinksMap->size() != remLinks->size()) {
			throw std::runtime_error(
					"Could not insert all removed links (are there duplicates?)");
		}

		// These two will be empty for now
		tpLinks = std::make_shared<EdgeContainerT>();
		tnLinks = std::make_shared<EdgeContainerT>();

		switch (posClass) {
		case TP:
			pos = tpLinks;
			break;
		case FN:
			pos = remLinks;
			break;
		case FP:
			pos = addLinks;
			break;
		case TN:
			pos = tnLinks;
			break;
		default:
			throw std::runtime_error("Unknown link class");
		}

		switch (negClass) {
		case TP:
			neg = tpLinks;
			break;
		case FN:
			neg = remLinks;
			break;
		case FP:
			neg = addLinks;
			break;
		case TN:
			neg = tnLinks;
			break;
		default:
			throw std::runtime_error("Unknown link class");
		}

		switch (posClass) {
		case TP:
			posStrmEndItSP = std::make_shared<TestEdgeIterator>(TPEIT,
					tpLinks->end(), eg->tpEnd(), eg->tnEnd());
			break;
		case FN:
			posStrmEndItSP = std::make_shared<TestEdgeIterator>(ECEIT,
					remLinks->end(), eg->tpEnd(), eg->tnEnd());
			break;
		case FP:
			posStrmEndItSP = std::make_shared<TestEdgeIterator>(ECEIT,
					addLinks->end(), eg->tpEnd(), eg->tnEnd());
			break;
		case TN:
			posStrmEndItSP = std::make_shared<TestEdgeIterator>(TNEIT,
					tnLinks->end(), eg->tpEnd(), eg->tnEnd());
			break;
		default:
			throw std::runtime_error("Unknown link class");
		}

		switch (negClass) {
		case TP:
			negStrmEndItSP = std::make_shared<TestEdgeIterator>(TPEIT,
					tpLinks->end(), eg->tpEnd(), eg->tnEnd());
			break;
		case FN:
			negStrmEndItSP = std::make_shared<TestEdgeIterator>(ECEIT,
					remLinks->end(), eg->tpEnd(), eg->tnEnd());
			break;
		case FP:
			negStrmEndItSP = std::make_shared<TestEdgeIterator>(ECEIT,
					addLinks->end(), eg->tpEnd(), eg->tnEnd());
			break;
		case TN:
			negStrmEndItSP = std::make_shared<TestEdgeIterator>(TNEIT,
					tnLinks->end(), eg->tpEnd(), eg->tnEnd());
			break;
		default:
			throw std::runtime_error("Unknown link class");
		}

		/**
		 eg->template generateTP(std::back_inserter(*tpLinks));
		 eg->template generateTN(std::back_inserter(*tnLinks));

		 tpGenerated = true;
		 tnGenerated = true;
		 */
	}

	/**
	 * Lock the dataset. No modifications can be done on the object after calling this method.
	 */
	void lock() {
		locked = true;
	}

	/**
	 * Generate positive instances.
	 */
	void genPos() {
		if (!locked) {
			switch (posClass) {
			case TP:
				if (!tpGenerated) {
					eg->template generateTP(std::back_inserter(*tpLinks));
					tpGenerated = true;
				}
				break;
			case TN:
				if (!tnGenerated) {
					eg->template generateTN(std::back_inserter(*tnLinks));
					tnGenerated = true;
				}
				break;
			default:
				break;
			}
		}
	}

	/**
	 * Generate negative instances.
	 */
	void genNeg() {
		if (!locked) {
			switch (negClass) {
			case TP:
				if (!tpGenerated) {
					eg->template generateTP(std::back_inserter(*tpLinks));
					tpGenerated = true;
				}
				break;
			case TN:
				if (!tnGenerated) {
					eg->template generateTN(std::back_inserter(*tnLinks));
					tnGenerated = true;
				}
				break;

			default:
				break;
			}
		}
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	TestData(TestData const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	TestData& operator =(TestData const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	TestData(TestData &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	TestData& operator =(TestData &&that) = default;

	/**
	 * @return The observed network.
	 */
	auto getObsNet() const {
		return obsNet;
	}

	/**
	 * @return The reference network.
	 */
	auto getRefNet() const {
		return refNet;
	}

	/**
	 * @return An iterator to the first positive link.
	 */
	auto posBegin() const {
		return pos->cbegin();
	}

	/**
	 * @return An iterator to one-past-the-last positive link.
	 */
	auto posEnd() const {
		return pos->cend();
	}

	/**
	 * @return An iterator to the first negative link.
	 */
	auto negBegin() const {
		return neg->cbegin();
	}

	/**
	 * @return An iterator to one-past-the-last negative link.
	 */
	auto negEnd() const {
		return neg->cend();
	}

	/**
	 * @return An iterator to the first positive link. The edges here may be streamed and not pre-generated.
	 */
	auto posStrmBegin() const {
		switch (posClass) {
		case TP:
			return TestEdgeIterator(TPEIT, tpLinks->end(), eg->tpBegin(),
					eg->tnEnd());
		case FN:
			return TestEdgeIterator(ECEIT, remLinks->begin(), eg->tpEnd(),
					eg->tnEnd());
		case FP:
			return TestEdgeIterator(ECEIT, addLinks->begin(), eg->tpEnd(),
					eg->tnEnd());
		case TN:
			return TestEdgeIterator(TNEIT, tpLinks->end(), eg->tpEnd(),
					eg->tnBegin());
		default:
			throw std::runtime_error("Unknown link class");
		}
	}

	/**
	 * @return An iterator to one-past-the-last positive link. The edges here may be streamed and not pre-generated.
	 */
	auto posStrmEnd() const {
		return *posStrmEndItSP;
	}

	/**
	 * @return An iterator to the first negative link. The edges here may be streamed and not pre-generated.
	 */
	auto negStrmBegin() const {
		switch (negClass) {
		case TP:
			return TestEdgeIterator(TPEIT, tpLinks->end(), eg->tpBegin(),
					eg->tnEnd());
		case FN:
			return TestEdgeIterator(ECEIT, remLinks->begin(), eg->tpEnd(),
					eg->tnEnd());
		case FP:
			return TestEdgeIterator(ECEIT, addLinks->begin(), eg->tpEnd(),
					eg->tnEnd());
		case TN:
			return TestEdgeIterator(TNEIT, tpLinks->end(), eg->tpEnd(),
					eg->tnBegin());
		default:
			throw std::runtime_error("Unknown link class");
		}
	}

	/**
	 * @return An iterator to one-past-the-last negative link. The edges here may be streamed and not pre-generated.
	 */
	auto negStrmEnd() const {
		return *negStrmEndItSP;
	}

	/**
	 * @return Number of positive links in the test set.
	 */
	std::size_t getNbPos() const {
		return pos->size();
	}

	/**
	 * @return Number of negative links in the test set.
	 */
	std::size_t getNbNeg() const {
		return neg->size();
	}

	/**
	 * @return The removed links in a set.
	 */
	std::shared_ptr<std::set<EdgeType> const> getRemLinksMap() const {
		return remLinksMap;
	}

	/**
	 * @return negClass.
	 */
	LinkClass getNegClass() const {
		return negClass;
	}

	/**
	 * @return posClass.
	 */
	LinkClass getPosClass() const {
		return posClass;
	}

	/**
	 * @return tnGenerated.
	 */
	bool isTnGenerated() const {
		return tnGenerated;
	}

	/**
	 * @return tpGenerated.
	 */
	bool isTpGenerated() const {
		return tpGenerated;
	}

	/**
	 * @return Whether the test data is locked.
	 */
	bool isLocked() const {
		return locked;
	}

	/**
	 * @return Edge generator.
	 */
	const std::shared_ptr<TestEdgeGenerator<NetworkT, EdgeContainerT> >& getEg() const {
		return eg;
	}

	/**
	 * Destructor.
	 */
	virtual ~TestData() = default;
};

/**
 * @brief Class to manipulate network by removing or adding edges.
 * @details
 * The following terminology is used:
 * When removing links:
 * TP: Edges of observed network
 * TN: Non-edges of reference network
 * FP: None
 * FN: Removed links
 * When adding links:
 * TP: Edges of reference network
 * TN: Non-edges of observed network
 * FP: Added links
 * FN: None
 * When adding/removing links:
 * TP: Edges of observed network - added links = Edges of reference network - removed links
 * TN: Non-edges of observed network - removed links = Non-edges of reference network - added links
 * FP: Added links
 * FN: Removed links
 * @tparam NetworkT The network type.
 */
template<typename NetworkT = UNetwork<>> class NetworkManipulator {

public:

	using NetworkSP = std::shared_ptr<NetworkT>; /**< Shared pointer to a network. */
	using NetworkCSP = std::shared_ptr<const NetworkT>; /**< Constant shared pointer to a network. */
	using NodeIdType = typename NetworkT::NodeIdType; /**< Nodes IDs type. */
	using LabelType = typename NetworkT::LabelType; /**< Nodes labels type. */
	using EdgeType = typename NetworkT::EdgeType; /**< The edge type. */
	template<typename ValueT> using NodeMap = typename NetworkT::template NodeMap<ValueT>; /**< Node map. */
	template<typename ValueT> using EdgeMap = typename NetworkT::template EdgeMap<ValueT>; /**< Edge map. */

protected:
	/**
	 * Creates test data by removing and/or adding edges to a network.
	 * The reference network is not modified. The two networks have the same external-internal ID mapping.
	 * @param refNet The reference network.
	 * @param obsNet The observed network.
	 * @param remLinks Removed links.
	 * @param addLinks Added links.
	 * @param remRatio Value between 0 and 1 that specifies the percentage of edges that are removed.
	 * @param addRatio Value between 0 and 1 that specifies the percentage of non-edges that are added.
	 * @param keepConnected Whether to keep the network connected.
	 * @param seed The random number generator's seed.
	 * @return The test data.
	 */
	static void createTestData(NetworkCSP refNet, NetworkSP obsNet,
			std::shared_ptr<std::vector<EdgeType>> remLinks,
			std::shared_ptr<std::vector<EdgeType>> addLinks, double remRatio,
			double addRatio, bool keepConnected, long int seed);

#ifdef WITH_OPENMP
	static bool parallel; /**< Enable/disable parallelism. */
#endif

public:

	/**
	 * Constructor.
	 */
	NetworkManipulator() = default;

	/**
	 * Generate a random spanning tree.
	 * @tparam InserterIterator Type of iterator to insert tree edges.
	 * @param net The network.
	 * @param seed The random number generator's seed.
	 * @param inserter An inserter iterator to insert the tree edges.
	 */
	template<typename InserterIterator> static void rst(NetworkCSP net,
			long int seed, InserterIterator inserter) {
		auto visited = net->template createNodeMap<uint8_t>(); // Unfortunately std::vector<bool> does not work as expected
		auto next = net->template createNodeMap<NodeIdType>();

		logger(logDebug, "Generating a random spanning tree...")
		RandomGen rng(seed);
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (auto it = visited.begin(); it < visited.end(); ++it) {
			*it = false;
		}

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (NodeIdType i = 0; i < net->getNbNodes(); i++) {
			next[i] = i;
		}

		visited[rng.getUInt(0, net->getNbNodes() - 1)] = true;
		for (NodeIdType i = 0; i < net->getNbNodes(); i++) {
			auto u = i;
			while (!visited[u]) {
				auto e = *Utilities::getRandom(net->neighborsBegin(u),
						net->neighborsEnd(u), rng.getInt());
				next[u] = NetworkT::end(e);
				u = next[u];
			}

			u = i;
			while (!visited[u]) {
				visited[u] = true;
				if (u != next[u]) {
					*inserter = NetworkT::makeEdge(u, next[u]);
					++inserter;
				}
				u = next[u];
			}
		}
		logger(logDebug, "Done")
	}

	/**
	 * Creates a new connected network from the current by extracting uniformly randomly a specified ratio of edges.
	 * The reference network is not modified. The two networks have the same external-internal ID mapping.
	 * @param net The reference network.
	 * @param ratio Value between 0 and 1 that specifies the percentage of edges that are removed.
	 * @param seed The random number generator's seed.
	 * @return A pair that contains a pointer to the resulting network and pointer to the set of extracted links.
	 */
	static std::pair<NetworkCSP, std::shared_ptr<std::vector<EdgeType>> > rndConExtract(
			NetworkCSP net, double ratio, long int seed);

	/**
	 * Creates a new network from the current by extracting uniformly randomly a specified ratio of edges.
	 * The reference network is not modified.  The two networks have the same external-internal ID mapping.
	 * @param net The reference network.
	 * @param ratio Value between 0 and 1 that specifies the percentage of edges that are removed.
	 * @param seed The random number generator's seed.
	 * @return A pair that contains a pointer to the resulting network and pointer to the set of extracted links.
	 */
	static std::pair<NetworkCSP, std::shared_ptr<std::vector<EdgeType>> > rndExtract(
			NetworkCSP net, double ratio, long int seed);

	/**
	 * Return all negative links except those passed in parameter. Be aware that this can be computationally intensive both in term of time and space.
	 * @param net The network.
	 * @param excepts The links that are excepted.
	 * @param inserter An inserter iterator.
	 */
	template<typename InserterIterator> static void getAllNegLinksExcept(
			NetworkCSP net, std::set<EdgeType> const &excepts,
			InserterIterator inserter) {
		for (auto it = net->nonEdgesBegin(); it != net->nonEdgesEnd(); ++it) {
			if ((excepts.count(*it) == 0)
					&& ((excepts.count(NetworkT::reverseEdge(*it)) == 0))) {
				*inserter = *it;
				++inserter;
			}
		}
	}

	/**
	 * Select random negative links except those passed in parameter.
	 * @param net The network.
	 * @param ratio The ratio of negative links to be selected.
	 * @param seed The random seed.
	 * @param excepts The links that are excepted.
	 * @param inserter The inserter iterator where the selected edges will be inserted.
	 */
	template<typename InserterIterator> static void getRndNegLinksExcept(
			NetworkCSP net, double ratio, long int seed,
			std::set<EdgeType> const &excepts, InserterIterator inserter) {
		Utilities::filter(net->rndNonEdgesBegin(ratio, seed),
				net->rndNonEdgesEnd(), excepts, inserter);
	}

	/**
	 * Select random negative links except those passed in parameter.
	 * @param net The network.
	 * @param ratio The ratio of positive links to be selected.
	 * @param seed The random seed.
	 * @param excepts The links that are excepted.
	 * @param inserter The inserter iterator where the selected edges will be inserted.
	 */
	template<typename InserterIterator> static void getRndPosLinksExcept(
			NetworkCSP net, double ratio, long int seed,
			std::set<EdgeType> const &excepts, InserterIterator inserter) {
		Utilities::filter(net->rndEdgesBegin(ratio, seed), net->rndEdgesEnd(),
				excepts, inserter);
	}

	// TODO combine the two methods below to offer the possibility of adding and removing edges at the same time.
	// This requires to change the container of remLink and addLinks to std::set, since a membership test
	// must be made to obtain true positive and true negative links.

	/**
	 * Creates test data by removing edges from a network.
	 * The reference network is not modified. The two networks have the same external-internal ID mapping.
	 * @param refNet The reference network.
	 * @param remRatio Value between 0 and 1 that specifies the percentage of edges that are removed.
	 * @param keepConnected Whether to keep the network connected.
	 * @param aTP Whether to use all true positive links in the test set.
	 * @param tpRatio Ratio of true positive inks to e used in the test set. This parameter is only relevant when aTP is false.
	 * @param aTN Whether to use all true negative links in the test set.
	 * @param tnRatio Ratio of true negative inks to e used in the test set. This parameter is only relevant when aTN is false.
	 * @param seed The random number generator's seed.
	 * @param preGenerateTPN Whether to pre-generate true positives and true negatives.
	 * @return The test data.
	 */
	static TestData<NetworkT, std::vector<EdgeType>> createTestDataRem(
			NetworkCSP refNet, double remRatio, bool keepConnected, bool aTP,
			double tpRatio, bool aTN, double tnRatio, long int seed,
			bool preGenerateTPN = true);

	/**
	 * Creates test data by adding edges to a network.
	 * The reference network is not modified. The two networks have the same external-internal ID mapping.
	 * @param refNet The reference network.
	 * @param addRatio Value between 0 and 1 that specifies the percentage of edges that are added.
	 * @param aTP Whether to use all true positive links in the test set.
	 * @param tpRatio Ratio of true positive inks to e used in the test set. This parameter is only relevant when aTP is false.
	 * @param aTN Whether to use all true negative links in the test set.
	 * @param tnRatio Ratio of true negative inks to e used in the test set. This parameter is only relevant when aTN is false.
	 * @param seed The random number generator's seed.
	 * @param preGenerateTPN Whether to pre-generate true positives and true negatives.
	 * @return The test data.
	 */
	static TestData<NetworkT, std::vector<EdgeType>> createTestDataAdd(
			NetworkCSP refNet, double addRatio, bool aTP, double tpRatio,
			bool aTN, double tnRatio, long int seed,
			bool preGenerateTPN = true);

	/**
	 * Creates test data by adding/removing edges from a network.
	 * The reference network is not modified. The two networks have the same external-internal ID mapping.
	 * @param refNet The reference network.
	 * @param remRatio Value between 0 and 1 that specifies the percentage of edges that are removed.
	 * @param addRatio Value between 0 and 1 that specifies the percentage of edges that are added.
	 * @param keepConnected Whether to keep the network connected.
	 * @param aTP Whether to use all true positive links in the test set.
	 * @param tpRatio Ratio of true positive inks to e used in the test set. This parameter is only relevant when aTP is false.
	 * @param aTN Whether to use all true negative links in the test set.
	 * @param tnRatio Ratio of true negative inks to e used in the test set. This parameter is only relevant when aTN is false.
	 * @param posClass Indicates which links will be considered the positive links.
	 * @param negClass Indicates which links will be considered the negative links.
	 * @param seed The random number generator's seed.
	 * @param preGenerateTPN Whether to pre-generate true positives and true negatives.
	 * @return The test data.
	 */
	static TestData<NetworkT, std::vector<EdgeType>> createTestData(
			NetworkCSP refNet, double remRatio, double addRatio,
			bool keepConnected, bool aTP, double tpRatio, bool aTN,
			double tnRatio, LinkClass posClass, LinkClass negClass,
			long int seed, bool preGenerateTPN = true);

	/**
	 * Creates test data from two networks. Only nodes common to both networks are considered.
	 * @param firstNet The first network.
	 * @param secondNet The second network.
	 * @param aTP Whether to use all true positive links in the test set.
	 * @param tpRatio Ratio of true positive inks to e used in the test set. This parameter is only relevant when aTP is false.
	 * @param aTN Whether to use all true negative links in the test set.
	 * @param tnRatio Ratio of true negative inks to e used in the test set. This parameter is only relevant when aTN is false.
	 * @param posClass Indicates which links will be considered the positive links.
	 * @param negClass Indicates which links will be considered the negative links.
	 * @param seed The random number generator's seed.
	 * @return The test data.
	 */
	static TestData<NetworkT, std::vector<EdgeType>> createTestDataSeqInter(
			NetworkCSP firstNet, NetworkCSP secondNet, bool aTP, double tpRatio,
			bool aTN, double tnRatio, LinkClass posClass, LinkClass negClass,
			long int seed);

	/**
	 * Creates test data from two networks. Nodes and edges not present in the second network are removed from the test set.
	 * @param firstNet The first network.
	 * @param secondNet The second network.
	 * @param aTP Whether to use all true positive links in the test set.
	 * @param tpRatio Ratio of true positive inks to e used in the test set. This parameter is only relevant when aTP is false.
	 * @param aTN Whether to use all true negative links in the test set.
	 * @param tnRatio Ratio of true negative inks to e used in the test set. This parameter is only relevant when aTN is false.
	 * @param posClass Indicates which links will be considered the positive links.
	 * @param negClass Indicates which links will be considered the negative links.
	 * @param seed The random number generator's seed.
	 * @param preGenerateTPN Whether to pre-generate true positives and true negatives.
	 * @return The test data.
	 */
	static TestData<NetworkT, std::vector<EdgeType>> createTestDataSeq(
			NetworkCSP firstNet, NetworkCSP secondNet, bool aTP, double tpRatio,
			bool aTN, double tnRatio, LinkClass posClass, LinkClass negClass,
			long int seed, bool preGenerateTPN = true);

	/**
	 * Check if a network is connected.
	 * @param net The network.
	 * @return True if the network net is connected, false otherwise.
	 */
	static bool isConnected(NetworkCSP net) {
		BFS<NetworkT> bfs(net);
		Collector<NetworkT> col;
		bfs.traverse(net->nodesBegin()->first, col);
		return col.getVisited().size() == net->getNbNodes();
	}

#ifdef WITH_OPENMP
	/**
	 * @return Whether parallelism is enabled.
	 */
	static bool isParallel() {
		return parallel;
	}

	/**
	 * Enable/disable parallelism.
	 * @param parallel True to enable parallelism, false to disable it.
	 */
	static void setParallel(bool parallel) {
		NetworkManipulator::parallel = parallel;
	}
#endif

	/**
	 * Destructor.
	 */
	virtual ~NetworkManipulator() = default;
};

}
/* namespace LinkPred */

#endif /* NETWORKMANIPULATOR_HPP_ */
