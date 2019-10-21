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
 * @brief Contains the implementation of a class to store and manage prediction results.
 */

#ifndef INCLUDE_PREDRESULTS_HPP_
#define INCLUDE_PREDRESULTS_HPP_

#include <linkpred/predictors/ulpredictor.hpp>
#include "linkpred/core/unetwork.hpp"
#include "linkpred/core/dnetwork.hpp"
#include "linkpred/perf/networkmanipulator.hpp"
#include "linkpred/utils/utilities.hpp"
#include <vector>
#include <memory>
#include <algorithm>
#include <functional>

namespace LinkPred {

/**
 * @brief A class to store and manage prediction results.
 * @tparam TestDataT The test data type.
 * @tparam LPredictorT The link predictor type.
 * @tparam ScoresContainerT The type of container storing scores.
 */
template<typename TestDataT = TestData<>, typename LPredictorT = ULPredictor<>,
		typename ScoresContainerT = std::vector<double>> class PredResults {

public:
	using ScoresIteratorT = typename ScoresContainerT::iterator; /**< Scores iterator type. */

protected:
	TestDataT testData; /**< Test data. */
	std::shared_ptr<LPredictorT> predictor; /**< The predictor. */
	ScoresContainerT posScores; /**< Scores of positive links. */
	ScoresContainerT negScores; /**< Scores of negative links. */
	bool posComputed = false; /**< Positive links scores computed? */
	bool negComputed = false; /**< Negative links scores computed? */
	bool topComputed = false; /**< Top links scores computed? */
	SortOrder posSortOrder = None; /**< Positive scores sorted? */
	SortOrder negSortOrder = None; /**< Negative scores sorted? */
	std::vector<typename TestDataT::EdgeType> topEdges; /**< To store top score edges. */
	ScoresContainerT topScores; /**< Scores of top links. */

public:

	/**
	 * @brief Test edges iterator.
	 */
	class ScoreIterator: public std::iterator<std::random_access_iterator_tag,
			const double, long int> {

		friend class PredResults; /**< PredResults is a friend. */

	protected:
		std::shared_ptr<LPredictorT> predictor; /**< The predictor. */
		typename TestDataT::TestEdgeIterator eit; /**< Edge iterator. */
		std::pair<typename TestDataT::EdgeType, double> esc; /**< Edge and score. */

	public:
		using pointer = typename std::iterator<std::random_access_iterator_tag, const std::pair<typename TestDataT::EdgeType, double>, long int>::pointer; /**< The pointer type associated with the iterator. */
		using reference = typename std::iterator<std::random_access_iterator_tag, const std::pair<typename TestDataT::EdgeType, double>, long int>::reference; /**< The reference type associated with the iterator. */
		using difference_type = typename std::iterator<std::random_access_iterator_tag, const std::pair<typename TestDataT::EdgeType, double>, long int>::difference_type; /**< The difference type associated with the iterator. */

		/**
		 * Constructor.
		 * @param td The test data object (owner).
		 */
		ScoreIterator(std::shared_ptr<LPredictorT> const & predictor,
				typename TestDataT::TestEdgeIterator const &eit) :
				predictor(predictor), eit(eit) {
		}

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		ScoreIterator(ScoreIterator const &that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		ScoreIterator& operator =(ScoreIterator const &that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		ScoreIterator(ScoreIterator &&that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		ScoreIterator& operator =(ScoreIterator &&that) = default;

		/**
		 * Dereference operator.
		 * @return A reference to the object to which the iterator points.
		 */
		reference operator*() {
			esc.first = *eit;
			esc.second = predictor->score(esc.first);
			return esc;
		}

		/**
		 * Arrow operator.
		 * @return A pointer to the object to which the iterator points.
		 */
		pointer operator->() {
			esc.first = *eit;
			esc.second = predictor->score(esc.first);
			return &esc;
		}

		/**
		 * Pre-increment operator.
		 * @return A reference to the new iterator.
		 */
		ScoreIterator& operator++() {
			++eit;
			return *this;
		}

		/**
		 * Pre-decrement operator.
		 * @return A reference to the new iterator.
		 */
		ScoreIterator& operator--() {
			--eit;
			return *this;
		}

		/**
		 * Post-increment operator.
		 * @return A reference to the new iterator.
		 */
		ScoreIterator operator++(int) {
			auto that = *this;
			++(*this);
			return that;
		}

		/**
		 * Post-decrement operator.
		 * @return A reference to the new iterator.
		 */
		ScoreIterator operator--(int) {
			auto that = *this;
			--(*this);
			return that;
		}

		/**
		 * Arithmetic + operator.
		 * @param n Increment value.
		 * @return The new iterator.
		 */
		ScoreIterator operator+(const difference_type &n) const {
			auto that = *this;
			that.eit += n;
			return that;
		}

		/**
		 * Arithmetic += operator.
		 * @param n Increment value.
		 * @return A reference to the new iterator.
		 */
		ScoreIterator& operator+=(const difference_type &n) {
			eit += n;
			return *this;
		}

		/**
		 * Arithmetic - operator.
		 * @param n Decrement value.
		 * @return The new iterator.
		 */
		ScoreIterator operator-(const difference_type &n) const {
			auto that = *this;
			that.eit -= n;
			return that;
		}

		/**
		 * Arithmetic -= operator.
		 * @param n Decrement value.
		 * @return A reference to the new iterator.
		 */
		ScoreIterator& operator-=(const difference_type &n) {
			eit -= n;
			return *this;
		}

		/**
		 * Difference between the present iterator and the one passed as parameter.
		 * @param that The other iterator.
		 * @return The difference between the current and that iterator.
		 */
		difference_type operator-(const ScoreIterator &that) const {
			return eit - that.eit;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this equals that.
		 */
		bool operator==(const ScoreIterator &that) const {
			return eit == that.eit;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is not equal to that.
		 */
		bool operator!=(const ScoreIterator &that) const {
			return !(*this == that);
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is less than that.
		 */
		bool operator<(const ScoreIterator &that) const {
			return eit < that.eit;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is greater that.
		 */
		bool operator>(const ScoreIterator &that) const {
			return eit > that.eit;
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is greater or equal to that.
		 */
		bool operator<=(const ScoreIterator &that) const {
			return !(*this > that);
		}

		/**
		 * @param that The other iterator.
		 * @return True if this is less or equal to that.
		 */
		bool operator>=(const ScoreIterator &that) const {
			return !(*this < that);
		}
	};

protected:
	std::shared_ptr<ScoreIterator> posSCoreStrmEndItSP; /**< End score iterator for negative links. Cached for performance. */
	std::shared_ptr<ScoreIterator> negSCoreStrmEndItSP; /**< End score iterator for negative links. Cached for performance. */

public:

	/**
	 * Constructor.
	 * @param testData The test data.
	 * @param predictor Link predictor.
	 */
	PredResults(TestDataT testData,
			std::shared_ptr<LPredictorT> const & predictor) :
			testData(testData), predictor(predictor) {
		if (!testData.isLocked()) {
			throw std::runtime_error("TestData must be locked");
		}

		if (testData.getEg()) {
			posSCoreStrmEndItSP = std::make_shared<ScoreIterator>(predictor,
					testData.posStrmEnd());
			negSCoreStrmEndItSP = std::make_shared<ScoreIterator>(predictor,
					testData.negStrmEnd());
		}
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	PredResults(PredResults const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	PredResults & operator =(PredResults const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	PredResults(PredResults && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	PredResults & operator =(PredResults && that) = default;

	/**
	 * Compute scores for positive links.
	 */
	void compPosScores() {
		if (!posComputed) {
//			testData.genPos(); // No need testData is locked
			posScores.resize(testData.getNbPos());
			predictor->predict(testData.posBegin(), testData.posEnd(),
					posScores.begin());
			Utilities::assertNoNaN(posScores.begin(), posScores.end());
			posComputed = true;
			posSortOrder = None;
		}
	}

	/**
	 * Compute scores for negative links.
	 */
	void compNegScores() {
		if (!negComputed) {
//			testData.genNeg();  // No need testData is locked
			negScores.resize(testData.getNbNeg());
			predictor->predict(testData.negBegin(), testData.negEnd(),
					negScores.begin());
			Utilities::assertNoNaN(negScores.begin(), negScores.end());
			negComputed = true;
			negSortOrder = None;
		}
	}

	/**
	 * Compute top scores. This method finds the negative edges with the highest scores.
	 * @param l Number of top links to compute.
	 */
	void compTopScores(std::size_t l) {
		if (l > testData.getObsNet()->getNbNonEdges()) {
			throw std::invalid_argument(
					"In TestData::compTopScores: l must be less or equal the number of negative links in the observed network.");
		}
		if (!topComputed) {
			topEdges.resize(l);
			topScores.resize(l);
			std::size_t k = predictor->top(l, topEdges.begin(),
					topScores.begin());
			topEdges.resize(k);
			topScores.resize(k);
			topComputed = true;
		}
	}

	/**
	 * @return Iterator to the first element in the scores of positive links.
	 */
	auto posBegin() const {
		return posScores.cbegin();
	}

	/**
	 * @return Iterator to one-past the last element in the scores of positive links.
	 */
	auto posEnd() const {
		return posScores.cend();
	}

	/**
	 * @return Iterator to the first element in the scores of negative links.
	 */
	auto negBegin() const {
		return negScores.cbegin();
	}

	/**
	 * @return Iterator to one-past the last element in the scores of negative links.
	 */
	auto negEnd() const {
		return negScores.cend();
	}

	/**
	 * @return Iterator to the first element in the scores of positive links. The scores here may be streamed and not pre-computed.
	 */
	auto posStrmBegin() const {
		return ScoreIterator(predictor, testData.posStrmBegin());
	}

	/**
	 * @return Iterator to one-past the last element in the scores of positive links. The scores here may be streamed and not pre-computed.
	 */
	auto posStrmEnd() const {
		return *posSCoreStrmEndItSP;
	}

	/**
	 * @return Iterator to the first element in the scores of negative links. The scores here may be streamed and not pre-computed.
	 */
	auto negStrmBegin() const {
		return ScoreIterator(predictor, testData.negStrmBegin());
	}

	/**
	 * @return Iterator to one-past the last element in the scores of negative links. The scores here may be streamed and not pre-computed.
	 */
	auto negStrmEnd() const {
		return *negSCoreStrmEndItSP;
	}

	/**
	 * @return Iterator to the first element in the top scores.
	 */
	auto topBegin() const {
		return topScores.cbegin();
	}

	/**
	 * @return Iterator to one-past the last element in the top scores.
	 */
	auto topEnd() const {
		return topScores.cend();
	}

	/**
	 * @return Iterator to the first element in the top edges.
	 */
	auto topEdgesBegin() const {
		return topEdges.cbegin();
	}

	/**
	 * @return Iterator to one-past the last element in the top edges.
	 */
	auto topEdgesEnd() const {
		return topEdges.cend();
	}

	/**
	 * @return Whether the negative scores have been computed.
	 */
	bool isNegComputed() const {
		return negComputed;
	}

	/**
	 * @return The sorting order of the negative scores.
	 */
	SortOrder getNegSortOrder() const {
		return negSortOrder;
	}

	/**
	 * Sort the negative scores.
	 * @param negSortOrder The request sorting order.
	 */
	void sortNeg(SortOrder negSortOrder) {
		if (this->negSortOrder != negSortOrder) {
			Utilities::sort(negScores.begin(), negScores.end(), negSortOrder);
			this->negSortOrder = negSortOrder;
		}
	}

	/**
	 * @return Whether the positive scores have been computed.
	 */
	bool isPosComputed() const {
		return posComputed;
	}

	/**
	 * @return The sorting order of the positive scores.		}

	 *
	 */
	SortOrder getPosSortOrder() const {
		return posSortOrder;
	}

	/**
	 * Sort the positive scores.
	 * @param posSortOrder The request sorting order.
	 */
	void sortPos(SortOrder posSortOrder) {
		if (this->posSortOrder != posSortOrder) {
			Utilities::sort(posScores.begin(), posScores.end(), posSortOrder);
			this->posSortOrder = posSortOrder;
		}
	}

	/**
	 * @return Whether the top scores have been computed.
	 */
	bool isTopComputed() const {
		return topComputed;
	}

	/**
	 * @return The link predictor.
	 */
	const std::shared_ptr<LPredictorT>& getPredictor() const {
		return predictor;
	}

	/**
	 * @return The test data.
	 */
	TestDataT & getTestData() {
		return testData;
	}

	/**
	 * @return The number of computed top edges.
	 */
	std::size_t getNbTop() const {
		if (!topComputed) {
			throw std::runtime_error("Top edges not yet computed");
		}
		return topEdges.size();
	}

	/**
	 * Destructor.
	 */
	virtual ~PredResults() = default;

};

} /* namespace LinkPred */

#endif /* INCLUDE_PREDRESULTS_HPP_ */
