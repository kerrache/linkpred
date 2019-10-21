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
 * @brief Contains the implementation of several performance measures.
 */

#ifndef INCLUDE_PERFMEASURE_HPP_
#define INCLUDE_PERFMEASURE_HPP_

#include "LinkPredConfig.hpp"
#include "linkpred/perf/predresults.hpp"
#include "linkpred/utils/log.hpp"
#include "linkpred/utils/utilities.hpp"
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <cmath>
#ifdef WITH_OPENMP
#include <omp.h>
#endif
#ifdef WITH_OPENMP
#include <omp.h>
#endif

namespace LinkPred {

using PerfResults = std::map<std::string, double>;
/**< Map of performance results.Performance results are stored in a map, with the key being the identifier of the performance measure and the data being the value. */

namespace PerfLambda {

using PerfLambdaT = std::function<double(std::size_t tp, std::size_t fn, std::size_t tn, std::size_t fp, std::size_t P, std::size_t N)>;
/**< Signature of performance lambdas. */

// Some useful lambdas.
/**
 * Recall.
 * @param tp The number of true positives.
 * @param fn The number of false negatives.
 * @param tn The number of true negatives.
 * @param fp The number of false positives.
 * @param P The number of positives.
 * @param N The number of negatives.
 */
auto rec = [](std::size_t tp, std::size_t fn, std::size_t tn, std::size_t fp,
		std::size_t P, std::size_t N) {
	return (double) tp / P;
};
/**
 * False positive rate.
 * @param tp The number of true positives.
 * @param fn The number of false negatives.
 * @param tn The number of true negatives.
 * @param fp The number of false positives.
 * @param P The number of positives.
 * @param N The number of negatives.
 */
auto fpr = [](std::size_t tp, std::size_t fn, std::size_t tn, std::size_t fp,
		std::size_t P, std::size_t N) {
	return (double) fp / N;
};
/**
 * Precision.
 * @param tp The number of true positives.
 * @param fn The number of false negatives.
 * @param tn The number of true negatives.
 * @param fp The number of false positives.
 * @param P The number of positives.
 * @param N The number of negatives.
 */
auto pre = [](std::size_t tp, std::size_t fn, std::size_t tn, std::size_t fp,
		std::size_t P, std::size_t N) {
	return (double) tp / (tp + fp);
};
/**
 * False negative rate.
 * @param tp The number of true positives.
 * @param fn The number of false negatives.
 * @param tn The number of true negatives.
 * @param fp The number of false positives.
 * @param P The number of positives.
 * @param N The number of negatives.
 */
auto fnr = [](std::size_t tp, std::size_t fn, std::size_t tn, std::size_t fp,
		std::size_t P, std::size_t N) {
	return (double) fn / P;
};
/**
 * True negative rate.
 * @param tp The number of true positives.
 * @param fn The number of false negatives.
 * @param tn The number of true negatives.
 * @param fp The number of false positives.
 * @param P The number of positives.
 * @param N The number of negatives.
 */
auto tnr = [](std::size_t tp, std::size_t fn, std::size_t tn, std::size_t fp,
		std::size_t P, std::size_t N) {
	return (double) tn / N;
};
/**
 * False omission rate.
 * @param tp The number of true positives.
 * @param fn The number of false negatives.
 * @param tn The number of true negatives.
 * @param fp The number of false positives.
 * @param P The number of positives.
 * @param N The number of negatives.
 */
auto fmr = [](std::size_t tp, std::size_t fn, std::size_t tn, std::size_t fp,
		std::size_t P, std::size_t N) {
	return (double) fn / (tn + fn);
};
/**
 * Accuracy.
 * @param tp The number of true positives.
 * @param fn The number of false negatives.
 * @param tn The number of true negatives.
 * @param fp The number of false positives.
 * @param P The number of positives.
 * @param N The number of negatives.
 */
auto acc = [](std::size_t tp, std::size_t fn, std::size_t tn, std::size_t fp,
		std::size_t P, std::size_t N) {
	return (double) (tp + tn) / (P + N);
};
/**
 * False discovery rate.
 * @param tp The number of true positives.
 * @param fn The number of false negatives.
 * @param tn The number of true negatives.
 * @param fp The number of false positives.
 * @param P The number of positives.
 * @param N The number of negatives.
 */
auto fdr = [](std::size_t tp, std::size_t fn, std::size_t tn, std::size_t fp,
		std::size_t P, std::size_t N) {
	return (double) fp / (tp + fp);
};
/**
 * Negative predictive value.
 * @param tp The number of true positives.
 * @param fn The number of false negatives.
 * @param tn The number of true negatives.
 * @param fp The number of false positives.
 * @param P The number of positives.
 * @param N The number of negatives.
 */
auto npv = [](std::size_t tp, std::size_t fn, std::size_t tn, std::size_t fp,
		std::size_t P, std::size_t N) {
	return (double) tn / (tn + fn);
};

}

/**
 * @brief Abstract performance measure.
 * @tparam PredResultsT The prediction results type.
 */
template<typename PredResultsT = PredResults<>> class PerfMeasure {

public:
	using ScoresIteratorT = typename PredResultsT::ScoresIteratorT; /**< Scores iterator type. */

protected:
	std::string name; /**< The name of the performance measure. */
#ifdef WITH_OPENMP
	bool parallel = false; /**< Enable/disable parallelism. */
#endif
#ifdef WITH_MPI
	bool distributed = false; /**< Enable/disable distributed parallelism. */
	MPI_Comm comm = MPI_COMM_WORLD; /**< The MPI communicator. */
#endif

public:
	/**
	 * Constructor.
	 */
	PerfMeasure() = default;

	/**
	 * Constructor.
	 * @param name The name of the performance measure.
	 */
	PerfMeasure(std::string const &name) :
			name(name) {

	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	PerfMeasure(PerfMeasure const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	PerfMeasure& operator =(PerfMeasure const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	PerfMeasure(PerfMeasure &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	PerfMeasure& operator =(PerfMeasure &&that) = default;

	/**
	 * @return The name of the performance measure.
	 */
	const std::string& getName() const {
		return name;
	}

	/**
	 * @param name The name of the performance measure.
	 */
	void setName(std::string const &name) {
		this->name = name;
	}

#ifdef WITH_OPENMP
	/**
	 * @return Whether parallelism is enabled.
	 */
	bool isParallel() const {
		return parallel;
	}

	/**
	 * Enable/disable parallelism.
	 * @param parallel True to enable parallelism, false to disable it.
	 */
	void setParallel(bool parallel) {
		this->parallel = parallel;
	}
#endif

#ifdef WITH_MPI
	/**
	 * @return Whether distributed memory parallelism is enabled.
	 */
	bool isDistributed() const {
		return distributed;
	}

	/**
	 * Enable/disable distributed memory parallelism.
	 * @param distributed True to enable distributed memory parallelism, false to disable it.
	 */
	void setDistributed(bool distributed) {
		this->distributed = distributed;
	}

	/**
	 * @return The MPI communicator.
	 */
	MPI_Comm getComm() const {
		return comm;
	}

	/**
	 * Set the MPI communicator.
	 * @param comm The new MPI communicator.
	 */
	void setComm(MPI_Comm const & comm) {
		this->comm = comm;
	}
#endif

	/**
	 * @return Whether the performance measure requires the generation of positive set. The default value is true.
	 */
	virtual bool requiresPos() const {
		return true;
	}

	/**
	 * @return Whether the performance measure requires the generation of negative set. The default value is true.
	 */
	virtual bool requiresNeg() const {
		return true;
	}

	/**
	 * @return Whether the performance measure requires network shuffling.
	 */
	virtual bool requiresShuffling() const {
		return false;
	}

	/**
	 * Computes the performance measure.
	 * @param predResults The prediction results.
	 * @param results To write results.
	 */
	virtual void eval(std::shared_ptr<PredResultsT> &predResults,
			PerfResults &results) = 0;

	/**
	 * Computes the performance measure.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param posSortOrder The sorting order of positive scores. This may be modified by the method.
	 * @param negSortOrder The sorting order of negative scores. This may be modified by the method.
	 * @param results To write results.
	 */
	virtual void eval(ScoresIteratorT posScoresBegin,
			ScoresIteratorT posScoresEnd, ScoresIteratorT negScoresBegin,
			ScoresIteratorT negScoresEnd, SortOrder &posSortOrder,
			SortOrder &negSortOrder, PerfResults &results) = 0;

	/**
	 * Destructor.
	 */
	virtual ~PerfMeasure() = default;

};

/**
 * @brief Abstract performance curve.
 * @tparam PredResultsT The prediction results type.
 */
template<typename PredResultsT = PredResults<>> class PerfCurve: public PerfMeasure<
		PredResultsT> {

public:
	using ScoresIteratorT = typename PerfMeasure<PredResultsT>::ScoresIteratorT; /**< Scores iterator type. */

protected:
	using PerfMeasure<PredResultsT>::name; /**< The name of the performance measure. */
#ifdef WITH_OPENMP
	using PerfMeasure<PredResultsT>::parallel; /**< Enable/disable parallelism. */
#endif
#ifdef WITH_MPI
	using PerfMeasure<PredResultsT>::comm; /**< The MPI communicator. */
	using PerfMeasure<PredResultsT>::distributed; /**< Enable/disable distributed parallelism. */
#endif

public:
	/**
	 * Constructor.
	 */
	PerfCurve() = default;

	/**
	 * Constructor.
	 * @param name The name of the performance measure.
	 */
	PerfCurve(std::string name) :
			PerfMeasure<PredResultsT>(name) {

	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	PerfCurve(PerfCurve const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	PerfCurve& operator =(PerfCurve const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	PerfCurve(PerfCurve &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	PerfCurve& operator =(PerfCurve &&that) = default;

	/**
	 * Computes the performance measure (typically the area under the curve or AUC).
	 * @param predResults The prediction results.
	 * @param results To write results.
	 */
	virtual void eval(std::shared_ptr<PredResultsT> &predResults,
			PerfResults &results) = 0;

	/**
	 * Computes the performance measure.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param posSortOrder The sorting order of positive scores. This may be modified by the method.
	 * @param negSortOrder The sorting order of negative scores. This may be modified by the method.
	 * @param results To write results.
	 */
	virtual void eval(ScoresIteratorT posScoresBegin,
			ScoresIteratorT posScoresEnd, ScoresIteratorT negScoresBegin,
			ScoresIteratorT negScoresEnd, SortOrder &posSortOrder,
			SortOrder &negSortOrder, PerfResults &results) = 0;

	/**
	 * Computes the performance curve.
	 * @param predResults The prediction results.
	 * @return A curve in the form of an std::vector of pairs representing the x and y coordinates.
	 */
	virtual std::vector<std::pair<double, double>> getCurve(
			std::shared_ptr<PredResultsT> &predResults) = 0;

	/**
	 * Computes the performance curve.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param posSortOrder The sorting order of positive scores. This may be modified by the method.
	 * @param negSortOrder The sorting order of negative scores. This may be modified by the method.
	 * @return A curve in the form of an std::vector of pairs representing the x and y coordinates.
	 */
	virtual std::vector<std::pair<double, double>> getCurve(
			ScoresIteratorT posScoresBegin, ScoresIteratorT posScoresEnd,
			ScoresIteratorT negScoresBegin, ScoresIteratorT negScoresEnd,
			SortOrder &posSortOrder, SortOrder &negSortOrder) = 0;

	/**
	 * Destructor.
	 */
	virtual ~PerfCurve() = default;

};

/**
 * @brief Receiver Operating  Characteristic curve.
 * @tparam PredResultsT The prediction results type.
 */
template<typename PredResultsT = PredResults<>> class ROC: public PerfCurve<
		PredResultsT> {
protected:
	bool strmEnabled = false; /**< Whether to enable streaming. */
	bool strmNeg = true; /**< Indicates which class to stream. If true negative scores are streamed, else positive scores are. Only used if enableStrm. */

public:
	using ScoresIteratorT = typename PerfMeasure<PredResultsT>::ScoresIteratorT; /**< Scores iterator type. */

protected:
	using PerfMeasure<PredResultsT>::name; /**< The name of the performance measure. */
#ifdef WITH_OPENMP
	using PerfMeasure<PredResultsT>::parallel; /**< Enable/disable parallelism. */
#endif
#ifdef WITH_MPI
	using PerfMeasure<PredResultsT>::comm; /**< The MPI communicator. */
	using PerfMeasure<PredResultsT>::distributed; /**< Enable/disable distributed parallelism. */
#endif

public:
	/**
	 * Constructor.
	 */
	ROC() {
		name = "ROC";
	}

	/**
	 * Constructor.
	 * @param name The name of the performance measure.
	 */
	ROC(std::string name) :
			PerfCurve<PredResultsT>(name) {

	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	ROC(ROC const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	ROC& operator =(ROC const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	ROC(ROC &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	ROC& operator =(ROC &&that) = default;

	/**
	 * @return Whether streaming is enabled.
	 */
	bool isStrmEnabled() const {
		return strmEnabled;
	}

	/**
	 * @param enableStrm Whether to enable streaming.
	 */
	void setStrmEnabled(bool strmEnabled) {
		this->strmEnabled = strmEnabled;
	}

	bool isStrmNeg() const {
		return strmNeg;
	}

	void setStrmNeg(bool strmNeg) {
		this->strmNeg = strmNeg;
	}

	/**
	 * Computes the area under the ROC curve.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param posSortOrder The sorting order of positive scores. This may be modified by the method.
	 * @param negSortOrder The sorting order of negative scores. This may be modified by the method.
	 * @param results To write results.
	 */
	virtual void eval(ScoresIteratorT posScoresBegin,
			ScoresIteratorT posScoresEnd, ScoresIteratorT negScoresBegin,
			ScoresIteratorT negScoresEnd, SortOrder &posSortOrder,
			SortOrder &negSortOrder, PerfResults &results) {
		logger(logDebug, "Computing ROC-AUC...")

		// Sort the smallest range as it leads to a better O()
		auto P = posScoresEnd - posScoresBegin;
		auto N = negScoresEnd - negScoresBegin;

		if (P <= N) {
			if (posSortOrder != Inc) {
				Utilities::sort(posScoresBegin, posScoresEnd, SortOrder::Inc);
				posSortOrder = Inc;
			}
		} else {
			if (negSortOrder != Inc) {
				Utilities::sort(negScoresBegin, negScoresEnd, SortOrder::Inc);
				negSortOrder = Inc;
			}
		}

#ifdef WITH_OPENMP
		results[name] = getROCAuc(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd, parallel);
#else
		results[name] = getROCAuc(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd);
#endif
		logger(logDebug, "Done")

	}

	/**
	 * Compute the area under the ROC curve. The smallest range must be sorted in increasing order (in case
	 * of a tie, the false negative ranges must be sorted).
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param parallel Whether to run in parallel.
	 * @return The area under the ROC curve.
	 */
	template<typename ScoresIterator> static double getROCAuc(
			ScoresIterator posScoresBegin, ScoresIterator posScoresEnd,
			ScoresIterator negScoresBegin, ScoresIterator negScoresEnd,
			bool parallel = false) {

		std::size_t nbLess = 0;
		std::size_t nbEq = 0;

		std::size_t P = posScoresEnd - posScoresBegin;
		std::size_t N = negScoresEnd - negScoresBegin;

		// The smallest range is sorted
		if (P <= N) {

#ifdef WITH_OPENMP
#pragma omp parallel for reduction(+ : nbLess), reduction(+ : nbEq) if (parallel)
#endif
			for (auto pit = negScoresBegin; pit < negScoresEnd; ++pit) {
				auto itl = std::lower_bound(posScoresBegin, posScoresEnd, *pit);
				auto itu = std::upper_bound(posScoresBegin, posScoresEnd, *pit);
				nbLess += posScoresEnd - itu;
				nbEq += itu - itl;
			}
		} else {

#ifdef WITH_OPENMP
#pragma omp parallel for reduction(+ : nbLess), reduction(+ : nbEq) if (parallel)
#endif
			for (auto pit = posScoresBegin; pit < posScoresEnd; ++pit) {
				auto itl = std::lower_bound(negScoresBegin, negScoresEnd, *pit);
				auto itu = std::upper_bound(negScoresBegin, negScoresEnd, *pit);
				nbLess += itl - negScoresBegin;
				nbEq += itu - itl;
			}

		}

		return (nbLess + 0.5 * nbEq) / (P * N);
	}

	/**
	 * Computes the area under the ROC curve.
	 * @param predResults The prediction results.
	 * @param results Iterator to write results.
	 */
	virtual void eval(std::shared_ptr<PredResultsT> &predResults,
			PerfResults &results) {

		logger(logDebug, "Computing ROC-AUC...")
		if (strmEnabled) {
			evalStream(predResults, results);
		} else {
			evalNoStream(predResults, results);
		}
		logger(logDebug, "Done")
	}

	/**
	 * Computes the area under the ROC curve.
	 * @param predResults The prediction results.
	 * @param results Iterator to write results.
	 */
	virtual void evalNoStream(std::shared_ptr<PredResultsT> &predResults,
			PerfResults &results) {

		logger(logDebug, "Computing ROC-AUC...")

		predResults->compNegScores();
		auto negScoresBegin = predResults->negBegin();
		auto negScoresEnd = predResults->negEnd();

		predResults->compPosScores();
		auto posScoresBegin = predResults->posBegin();
		auto posScoresEnd = predResults->posEnd();

		// Sort the smallest range as it leads to a better O()
		auto P = posScoresEnd - posScoresBegin;
		auto N = negScoresEnd - negScoresBegin;

		if (P <= N) {
			predResults->sortPos(SortOrder::Inc);
		} else {
			predResults->sortNeg(SortOrder::Inc);
		}

#ifdef WITH_OPENMP
		results[name] = getROCAuc(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd, parallel);
#else
		results[name] = getROCAuc(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd);
#endif
		logger(logDebug, "Done")
	}

	/**
	 * Compute the area under the ROC curve with streaming. The smallest range must be sorted in increasing order (in case
	 * of a tie, the false negative ranges must be sorted).
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param parallel Whether to run in parallel.
	 * @return The area under the ROC curve.
	 */
	template<typename StrmScoresIterator, typename StoredScoresIterator> double getROCAucStrm(
			StrmScoresIterator strmBegin, StrmScoresIterator strmEnd,
			StoredScoresIterator strdBegin, StoredScoresIterator strdEnd) {

		std::size_t nbStrm = 0;
		auto beginSP = std::make_shared<StrmScoresIterator>(strmBegin);
		auto endSP = std::make_shared<StrmScoresIterator>(strmEnd);
#ifdef WITH_MPI
		int nbProcs = 1;
		int procID = 0;
		if (distributed) {
			MPI_Comm_size(comm, &nbProcs);
			MPI_Comm_rank(comm, &procID);
			nbStrm = strmEnd - strmBegin;
			auto lr = Utilities::localRange(nbStrm, nbProcs, procID);
			beginSP = std::make_shared<StrmScoresIterator>(strmBegin + lr.first);
			endSP = std::make_shared<StrmScoresIterator>(strmBegin + lr.second);
		}
#endif

		auto begin = *beginSP;
		auto end = *endSP;

		std::size_t nbLess = 0;
		std::size_t nbEq = 0;

		std::size_t nbStrd = strdEnd - strdBegin;

		if (nbStrm == 0) {
#ifdef WITH_OPENMP
#pragma omp parallel for reduction(+ : nbStrm), reduction(+ : nbLess), reduction(+ : nbEq) if (parallel)
		for (auto pit = begin; pit < end; ++pit) {
#else
			for (auto pit = begin; pit != end; ++pit) {
#endif
				double sc = pit->second;
				auto itl = std::lower_bound(strdBegin, strdEnd, sc);
				auto itu = std::upper_bound(strdBegin, strdEnd, sc);
				nbStrm++;
				nbLess += strdEnd - itu;
				nbEq += itu - itl;
			}
		} else {
#ifdef WITH_OPENMP
#pragma omp parallel for reduction(+ : nbLess), reduction(+ : nbEq) if (parallel)
		for (auto pit = begin; pit < end; ++pit) {
#else
			for (auto pit = begin; pit != end; ++pit) {
#endif
				double sc = pit->second;
				auto itl = std::lower_bound(strdBegin, strdEnd, sc);
				auto itu = std::upper_bound(strdBegin, strdEnd, sc);
				nbLess += strdEnd - itu;
				nbEq += itu - itl;
			}
		}

#ifdef WITH_MPI
		if (distributed) {
			unsigned long long lNbLessEq[2];
			unsigned long long gNbLessEq[2];
			lNbLessEq[0] = nbLess;
			lNbLessEq[1] = nbEq;
			MPI_Reduce(lNbLessEq, gNbLessEq, 2, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
					0, comm);
			nbLess = gNbLessEq[0];
			nbEq = gNbLessEq[1];
		}
#endif
		return (nbLess + 0.5 * nbEq) / (nbStrd * nbStrm);
	}

	/**
	 * Computes the area under the ROC curve with streaming scores.
	 * @param predResults The prediction results.
	 * @param results Iterator to write results.
	 */
	virtual void evalStream(std::shared_ptr<PredResultsT> &predResults,
			PerfResults &results) {
		if (strmNeg) {
			predResults->compPosScores();
			auto posScoresBegin = predResults->posBegin();
			auto posScoresEnd = predResults->posEnd();
			predResults->sortPos(SortOrder::Inc);
			results[name] = getROCAucStrm(predResults->negStrmBegin(),
					predResults->negStrmEnd(), posScoresBegin, posScoresEnd);
		} else {
			predResults->compNegScores();
			auto negScoresBegin = predResults->negBegin();
			auto negScoresEnd = predResults->negEnd();
			predResults->sortNeg(SortOrder::Inc);
			results[name] = 1.0
					- getROCAucStrm(predResults->posStrmBegin(),
							predResults->posStrmEnd(), negScoresBegin,
							negScoresEnd);
		}
	}

	/**
	 * Compute the threshold of the ROC curve. Ranges must be sorted in increasing order.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 */
	template<typename ScoresIterator> static std::vector<double> getThresholds(
			ScoresIterator posScoresBegin, ScoresIterator posScoresEnd,
			ScoresIterator negScoresBegin, ScoresIterator negScoresEnd) {

		std::size_t P = posScoresEnd - posScoresBegin;
		std::size_t N = negScoresEnd - negScoresBegin;

		std::vector<double> allScores;
		allScores.resize(P + N);
		std::merge(posScoresBegin, posScoresEnd, negScoresBegin, negScoresEnd,
				allScores.begin());

		std::vector<double> unique;
		for (auto it = allScores.cbegin(); it < allScores.cend();) {
			unique.push_back(*it);
			auto rit = it + 1;
			while ((rit < posScoresEnd) && (*rit == *it)) {
				++rit;
			}
			it = rit;
		}
		return unique;
	}

	/**
	 * Compute the ROC curve. Both ranges must be sorted.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param parallel Whether to run in parallel.
	 */
	template<typename ScoresIterator> static std::vector<
			std::pair<double, double>> getROCCurve(
			ScoresIterator posScoresBegin, ScoresIterator posScoresEnd,
			ScoresIterator negScoresBegin, ScoresIterator negScoresEnd,
			bool parallel = false) {

		std::size_t P = posScoresEnd - posScoresBegin;
		std::size_t N = negScoresEnd - negScoresBegin;

		auto threshs = getThresholds(posScoresBegin, posScoresEnd,
				negScoresBegin, negScoresEnd);

		std::vector<std::pair<double, double>> curve;
		curve.resize(threshs.size() + 1);
		std::size_t maxInd = threshs.size();

//		std::cout << "TP\tFP\n";
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < threshs.size(); i++) {

			auto itpl = std::lower_bound(posScoresBegin, posScoresEnd,
					threshs[i]);
			std::size_t tp = posScoresEnd - itpl;

			auto itnl = std::lower_bound(negScoresBegin, negScoresEnd,
					threshs[i]);
			std::size_t fp = negScoresEnd - itnl;

//			std::cout << tp << "\t" << fp << std::endl;

			curve[maxInd - i].first = (double) fp / N;
			curve[maxInd - i].second = (double) tp / P;
		}
//		std::cout << "-----------\n";

// First entry (boundary)
		curve[0].first = 0;
		curve[0].second = 0;
		return curve;
	}

	/**
	 * Computes the ROC curve.
	 * @param predResults The prediction results.
	 * @return A curve in the form of an std::vector of pairs representing the x and y coordinates.
	 */
	virtual std::vector<std::pair<double, double>> getCurve(
			std::shared_ptr<PredResultsT> &predResults) {

		predResults->compNegScores();
		predResults->sortNeg(SortOrder::Inc);
		auto negScoresBegin = predResults->negBegin();
		auto negScoresEnd = predResults->negEnd();

		predResults->compPosScores();
		predResults->sortPos(SortOrder::Inc);
		auto posScoresBegin = predResults->posBegin();
		auto posScoresEnd = predResults->posEnd();

#ifdef WITH_OPENMP
		return getROCCurve(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd, parallel);
#else
		return getROCCurve(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd);
#endif
	}

	/**
	 * Computes the ROC curve.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param posSortOrder The sorting order of positive scores. This may be modified by the method.
	 * @param negSortOrder The sorting order of negative scores. This may be modified by the method.
	 * @return A curve in the form of an std::vector of pairs representing the x and y coordinates.
	 */
	virtual std::vector<std::pair<double, double>> getCurve(
			ScoresIteratorT posScoresBegin, ScoresIteratorT posScoresEnd,
			ScoresIteratorT negScoresBegin, ScoresIteratorT negScoresEnd,
			SortOrder &posSortOrder, SortOrder &negSortOrder) {

		if (posSortOrder != Inc) {
			Utilities::sort(posScoresBegin, posScoresEnd, SortOrder::Inc);
			posSortOrder = Inc;
		}
		if (negSortOrder != Inc) {
			Utilities::sort(negScoresBegin, negScoresEnd, SortOrder::Inc);
			negSortOrder = Inc;
		}

#ifdef WITH_OPENMP
		return getROCCurve(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd, parallel);
#else
		return getROCCurve(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd);
#endif
	}

	/**
	 * Destructor.
	 */
	virtual ~ROC() = default;
};

/**
 * @brief The precision recall curve.
 * @tparam PredResultsT The prediction results type.
 */
template<typename PredResultsT = PredResults<>> class PR: public PerfCurve<
		PredResultsT> {

public:
	using ScoresIteratorT = typename PerfMeasure<PredResultsT>::ScoresIteratorT; /**< Scores iterator type. */

protected:
	using PerfMeasure<PredResultsT>::name; /**< The name of the performance measure. */
#ifdef WITH_OPENMP
	using PerfMeasure<PredResultsT>::parallel; /**< Enable/disable parallelism. */
#endif

public:

	/**
	 * @brief Interpolation methods used for the computation of the PR-AUC.
	 */
	enum InterpolMethod {
		LIN, /**< Linear interpolation (Trapezoidal rule). */
		DGI /**< Davis-Goadrich nonlinear interpolation. */
	};

protected:
	InterpolMethod interpolMethod = DGI; /**< Interpolation method used for the computation of the PR-AUC. */

public:

	/**
	 * Constructor.
	 */
	PR() {
		name = "PR";
	}

	/**
	 * Constructor.
	 * @param name The name of the performance measure.
	 */
	PR(std::string name) :
			PerfCurve<PredResultsT>(name) {

	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	PR(PR const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	PR& operator =(PR const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	PR(PR &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	PR& operator =(PR &&that) = default;

	/**
	 * Compute the area under the Precision-Recall Curve using linear interpolation (trapezoidal rule).
	 * Notice that this method may over-estimate the actual area especially when the points on the
	 * curve are distant from one-another.
	 * @param tpsBegin Iterator to the first true positive count.
	 * @param tpsEnd Iterator to one-past-the-last true positive count.
	 * @param fpsBegin Iterator to the first false positive count.
	 * @param fpsEnd Iterator to one-past-the-last false positive count.
	 * @param P The number of positive instances.
	 * @param bcz Whether the boundary value of precision is zero.
	 * @param parallel Whether to run in parallel.
	 */
	template<typename CountIterator> static double getPRAucLIN(
			CountIterator tpsBegin, CountIterator tpsEnd,
			CountIterator fpsBegin, CountIterator fpsEnd, std::size_t P,
			bool bcz, bool parallel = false) {

		std::size_t nbtps = tpsEnd - fpsBegin;
		std::size_t nbfps = fpsEnd - fpsBegin;
		if (nbtps != nbfps) {
			throw std::runtime_error(
					"True positive and false positive count ranges must be equal");
		}
		if (nbtps == 0) {
			throw std::runtime_error("Ranges must not be empty");
		}

		double prAucLIN = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for reduction(+ : prAucLIN) if (parallel)
#endif
		for (std::size_t i = 0; i < nbtps - 1; i++) {
			double tpsi = *(tpsBegin + i);
			double tpsj = *(tpsBegin + i + 1);
			double fpsi = *(fpsBegin + i);
			double fpsj = *(fpsBegin + i + 1);

			prAucLIN += 0.5 * (tpsi - tpsj)
					* (tpsj / (tpsj + fpsj) + tpsi / (tpsi + fpsi));
		}

		// Last entry
		std::size_t last = nbtps - 1;
		auto tpsl = *(tpsBegin + last);
		auto fpsl = *(fpsBegin + last);
		if (bcz) {
			prAucLIN += 0.5 * tpsl * (tpsl / (tpsl + fpsl)); // Here: pre[pre.size() - 1] = 0;
		} else {
			prAucLIN += tpsl * (tpsl / (tpsl + fpsl)); // Here: pre[pre.size() - 1] = pre[pre.size() - 2];
		}

		return prAucLIN / P;
	}

	/**
	 * Compute the area under the Precision-Recall Curve using Davis-Goadrich nonlinear interpolation.
	 * See: Jesse Davis and Mark Goadrich (2006) The relationship between Precision-Recall and ROC curves.
	 * Proceedings of the 23rd international conference on Machine learning. pp.233–240
	 * @param tpsBegin Iterator to the first true positive count.
	 * @param tpsEnd Iterator to one-past-the-last true positive count.
	 * @param fpsBegin Iterator to the first false positive count.
	 * @param fpsEnd Iterator to one-past-the-last false positive count.
	 * @param P The number of positive instances.
	 * @param bcz Whether the boundary value of precision is zero.
	 * @param parallel Whether to run in parallel.
	 */
	template<typename CountIterator> static double getPRAucDGI(
			CountIterator tpsBegin, CountIterator tpsEnd,
			CountIterator fpsBegin, CountIterator fpsEnd, std::size_t P,
			bool bcz, bool parallel = false) {

		std::size_t nbtps = tpsEnd - tpsBegin;
		std::size_t nbfps = tpsEnd - tpsBegin;
		if (nbtps != nbfps) {
			throw std::runtime_error(
					"True positive and false positive count ranges must be equal");
		}
		if (nbtps == 0) {
			throw std::runtime_error("Ranges must not be empty");
		}

		double prAucDGI = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for reduction(+ : prAucDGI) if (parallel)
#endif
		for (std::size_t i = 0; i < nbtps - 1; i++) {
			auto tpsi = *(tpsBegin + i);
			auto tpsj = *(tpsBegin + i + 1);
			auto fpsi = *(fpsBegin + i);
			auto fpsj = *(fpsBegin + i + 1);

			if (tpsi != tpsj) {
				double a = tpsj;
				double b = 1 + (double) (fpsi - fpsj) / (tpsi - tpsj);
				double c = tpsj + fpsj;
				double d = tpsi - tpsj;

				double v = -(c * std::log(std::abs(b * d + c))) / (b * b)
						+ (a * std::log(std::abs(b * d + c))) / b + d / b
						+ (c * std::log(std::abs(c))) / (b * b)
						- (a * std::log(std::abs(c))) / b;
				prAucDGI += v;
			}
		}
		// Last entry
		std::size_t last = nbtps - 1;
		double tpsl = *(tpsBegin + last);
		double fpsl = *(fpsBegin + last);

		if (bcz) {
			double b = 1 + fpsl / tpsl;
			double d = tpsl;
			prAucDGI += d / b; // Here: pre[pre.size() - 1] = 0;
		} else {
			prAucDGI += tpsl * (tpsl / (tpsl + fpsl)); // Here: pre[pre.size() - 1] = pre[pre.size() - 2];
		}

		prAucDGI /= P;
		return prAucDGI;
	}

	/**
	 * Compute the area under the Precision-Recall Curve. Ranges must be sorted in increasing order.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param interpolMethod Interpolation method used in the integral.
	 * @param parallel Whether to run in parallel.
	 */
	template<typename ScoresIterator> static double getPRAuc(
			ScoresIterator posScoresBegin, ScoresIterator posScoresEnd,
			ScoresIterator negScoresBegin, ScoresIterator negScoresEnd,
			InterpolMethod interpolMethod, bool parallel = false) {

//		std::cerr << posScoresEnd - posScoresBegin << " ---- "
//				<< negScoresEnd - negScoresBegin << std::endl;

		std::size_t P = posScoresEnd - posScoresBegin;

		auto threshs = getThresholds(posScoresBegin, posScoresEnd,
				negScoresBegin, negScoresEnd);

		std::vector<std::size_t> tps;
		tps.resize(threshs.size());
		std::vector<std::size_t> fps;
		fps.resize(threshs.size());

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < threshs.size(); i++) {

			auto itpl = std::lower_bound(posScoresBegin, posScoresEnd,
					threshs[i]);
			tps[i] = posScoresEnd - itpl;

			auto itnl = std::lower_bound(negScoresBegin, negScoresEnd,
					threshs[i]);
			fps[i] = negScoresEnd - itnl;
		}

		// Last entry (boundary)
		bool bcz;
		std::size_t last = threshs.size() - 1;
		auto itnu = std::upper_bound(negScoresBegin, negScoresEnd,
				threshs[last]);
		if (itnu != negScoresEnd) {
			bcz = true; // Here: pre[pre.size() - 1] = 0;
		} else {
			bcz = false; // Here: pre[pre.size() - 1] = pre[pre.size() - 2];
		}

		switch (interpolMethod) {

		case LIN:
			return getPRAucLIN(tps.cbegin(), tps.cend(), fps.cbegin(),
					fps.cend(), P, bcz, parallel);

		case DGI:
			return getPRAucDGI(tps.cbegin(), tps.cend(), fps.cbegin(),
					fps.cend(), P, bcz, parallel);

		default:
			throw std::runtime_error("Unknown interpolation method");
		}
	}

	/**
	 * Computes the performance measure.
	 * @param predResults The prediction results.
	 * @param results Iterator to write results.
	 */
	virtual void eval(std::shared_ptr<PredResultsT> &predResults,
			PerfResults &results) {

		logger(logDebug, "Computing PR-AUC...")

		predResults->compNegScores();
		predResults->sortNeg(SortOrder::Inc);
		auto negScoresBegin = predResults->negBegin();
		auto negScoresEnd = predResults->negEnd();

		predResults->compPosScores();
		predResults->sortPos(SortOrder::Inc);
		auto posScoresBegin = predResults->posBegin();
		auto posScoresEnd = predResults->posEnd();

#ifdef WITH_OPENMP
		results[name] = getPRAuc(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd, interpolMethod, parallel);
#else
		results[name] = getPRAuc(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd, interpolMethod);
#endif

		logger(logDebug, "Done")
	}

	/**
	 * Computes the area under the PR curve.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param posSortOrder The sorting order of positive scores. This may be modified by the method.
	 * @param negSortOrder The sorting order of negative scores. This may be modified by the method.
	 * @param results To write results.
	 */
	virtual void eval(ScoresIteratorT posScoresBegin,
			ScoresIteratorT posScoresEnd, ScoresIteratorT negScoresBegin,
			ScoresIteratorT negScoresEnd, SortOrder &posSortOrder,
			SortOrder &negSortOrder, PerfResults &results) {
		logger(logDebug, "Computing PR-AUC...")

		if (posSortOrder != Inc) {
			Utilities::sort(posScoresBegin, posScoresEnd, SortOrder::Inc);
			posSortOrder = Inc;
		}
		if (negSortOrder != Inc) {
			Utilities::sort(negScoresBegin, negScoresEnd, SortOrder::Inc);
			negSortOrder = Inc;
		}

#ifdef WITH_OPENMP
		results[name] = getPRAuc(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd, interpolMethod, parallel);
#else
		results[name] = getPRAuc(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd, interpolMethod);
#endif
		logger(logDebug, "Done")

	}

	/**
	 * Compute the threshold of the Precision-Recall Curve. Ranges must be sorted in increasing order.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 */
	template<typename ScoresIterator> static std::vector<double> getThresholds(
			ScoresIterator posScoresBegin, ScoresIterator posScoresEnd,
			ScoresIterator negScoresBegin, ScoresIterator negScoresEnd) {

		std::size_t P = posScoresEnd - posScoresBegin;
		std::size_t N = negScoresEnd - negScoresBegin;

		std::vector<double> allScores;
		allScores.resize(P + N);
		std::merge(posScoresBegin, posScoresEnd, negScoresBegin, negScoresEnd,
				allScores.begin());

		std::vector<double> unique;
		for (auto it = allScores.cbegin(); it < allScores.cend();) {
			unique.push_back(*it);
			auto rit = it + 1;
			while ((rit < posScoresEnd) && (*rit == *it)) {
				++rit;
			}
			it = rit;
		}
		return unique;
	}

	/**
	 * Compute the Precision-Recall Curve. Ranges must sorted in increasing order.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param parallel Whether to run in parallel.
	 */
	template<typename ScoresIterator> static std::vector<
			std::pair<double, double>> getPRCurve(ScoresIterator posScoresBegin,
			ScoresIterator posScoresEnd, ScoresIterator negScoresBegin,
			ScoresIterator negScoresEnd, bool parallel = false) {

		std::size_t P = posScoresEnd - posScoresBegin;

		auto threshs = getThresholds(posScoresBegin, posScoresEnd,
				negScoresBegin, negScoresEnd);

		std::vector<std::pair<double, double>> curve;
		curve.resize(threshs.size() + 1);
		std::size_t maxInd = threshs.size();

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < threshs.size(); i++) {

			auto itpl = std::lower_bound(posScoresBegin, posScoresEnd,
					threshs[i]);
			std::size_t tp = posScoresEnd - itpl;

			auto itnl = std::lower_bound(negScoresBegin, negScoresEnd,
					threshs[i]);
			std::size_t fp = negScoresEnd - itnl;

			curve[maxInd - i].first = (double) tp / P;
			curve[maxInd - i].second = (double) tp / (tp + fp);
		}

		// First entry (boundary)
		curve[0].first = 0;
		std::size_t last = threshs.size() - 1;
		auto itnu = std::upper_bound(negScoresBegin, negScoresEnd,
				threshs[last]);
		if (itnu != negScoresEnd) {
			curve[0].second = 0; // Here: pre[pre.size() - 1] = 0;
		} else {
			curve[0].second = curve[1].second; // Here: pre[pre.size() - 1] = pre[pre.size() - 2];
		}

		return curve;
	}

	/**
	 * Compute the PR curve.
	 * @param predResults The prediction results.
	 * @return A curve in the form of an std::vector of pairs representing the x and y coordinates.
	 */
	virtual std::vector<std::pair<double, double>> getCurve(
			std::shared_ptr<PredResultsT> &predResults) {

		predResults->compNegScores();
		predResults->sortNeg(SortOrder::Inc);
		auto negScoresBegin = predResults->negBegin();
		auto negScoresEnd = predResults->negEnd();

		predResults->compPosScores();
		predResults->sortPos(SortOrder::Inc);
		auto posScoresBegin = predResults->posBegin();
		auto posScoresEnd = predResults->posEnd();

#ifdef WITH_OPENMP
		return getPRCurve(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd, parallel);
#else
		return getPRCurve(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd);
#endif
	}

	/**
	 * Compute the PR curve.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param posSortOrder The sorting order of positive scores. This may be modified by the method.
	 * @param negSortOrder The sorting order of negative scores. This may be modified by the method.
	 * @return A curve in the form of an std::vector of pairs representing the x and y coordinates.
	 */
	virtual std::vector<std::pair<double, double>> getCurve(
			ScoresIteratorT posScoresBegin, ScoresIteratorT posScoresEnd,
			ScoresIteratorT negScoresBegin, ScoresIteratorT negScoresEnd,
			SortOrder &posSortOrder, SortOrder &negSortOrder) {

		if (posSortOrder != Inc) {
			Utilities::sort(posScoresBegin, posScoresEnd, SortOrder::Inc);
			posSortOrder = Inc;
		}
		if (negSortOrder != Inc) {
			Utilities::sort(negScoresBegin, negScoresEnd, SortOrder::Inc);
			negSortOrder = Inc;
		}

#ifdef WITH_OPENMP
		return getPRCurve(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd, parallel);
#else
		return getPRCurve(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd);
#endif
	}

	/**
	 * @return The interpolation method.
	 */
	InterpolMethod getInterpolMethod() const {
		return interpolMethod;
	}

	/**
	 * Set the interpolation method.
	 * @param interpolMethod The new value of the interpolation method.
	 */
	void setInterpolMethod(InterpolMethod interpolMethod) {
		this->interpolMethod = interpolMethod;
	}

	/**
	 * Destructor.
	 */
	virtual ~PR() = default;

};

/**
 * @brief Compute top precision.
 * @tparam PredResultsT The prediction results type.
 */
template<typename PredResultsT = PredResults<>> class TPR: public PerfMeasure<
		PredResultsT> {

public:
	using ScoresIteratorT = typename PerfMeasure<PredResultsT>::ScoresIteratorT; /**< Scores iterator type. */

protected:
	using PerfMeasure<PredResultsT>::name; /**< The name of the performance measure. */
#ifdef WITH_OPENMP
	using PerfMeasure<PredResultsT>::parallel; /**< Enable/disable parallelism. */
#endif

protected:
	std::size_t l; /**< The number of links to consider. */
	bool useTopMethod = false; /**< Whether to use the top method of the link predictor. */

public:
	/**
	 * Constructor.
	 * @param l The number of links to consider.
	 */
	TPR(std::size_t l) :
			l(l) {
		name = "TPR";
	}

	/**
	 * Constructor.
	 * @param l The number of links to consider.
	 * @param name The name of the performance measure.
	 */
	TPR(std::size_t l, std::string name) :
			PerfMeasure<PredResultsT>(name), l(l) {
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	TPR(TPR const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	TPR& operator =(TPR const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	TPR(TPR &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	TPR& operator =(TPR &&that) = default;

	/**
	 * Ranges must be sorted in decreasing order.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param l The number of links to consider.
	 */
	template<typename ScoresIterator> static double getTPR(
			ScoresIterator posScoresBegin, ScoresIterator posScoresEnd,
			ScoresIterator negScoresBegin, ScoresIterator negScoresEnd,
			std::size_t l) {

		if (l
				> static_cast<std::size_t>(std::distance(posScoresBegin,
						posScoresEnd)
						+ std::distance(negScoresBegin, negScoresEnd))) {
			throw std::invalid_argument(
					"l must be less or equal the total number of testset elements");
		}
		if (l
				> static_cast<std::size_t>(std::distance(posScoresBegin,
						posScoresEnd))) {
			std::cerr
					<< "# Warning in getTPR: l is larger than the number of positive links. The TPR cannot reach 1 in this case!\n";
		}

		double pj = 0;
		auto cpit = posScoresEnd - 1;
		auto cnit = negScoresEnd - 1;
		std::size_t sel = 0;
		while (sel < l) {
			double ps = *cpit;
			double ns = *cnit;
			if (ps > ns) {
				pj++;
				sel++;
				cpit--;
			} else if (ps < ns) {
				sel++;
				cnit--;
			} else {
				auto plit = std::lower_bound(posScoresBegin, cpit, ps);
				auto nlit = std::lower_bound(negScoresBegin, cnit, ns);
				std::size_t nps = cpit - plit + 1;
				std::size_t nns = cnit - nlit + 1;
				std::size_t rem = l - sel;
				if (rem <= nps + nns) {
					pj += (double) rem * nps / (nps + nns);
					sel = l;
				} else {
					pj += nps;
					sel += nps + nns;
					cpit = plit - 1;
					cnit = nlit - 1;
				}
			}
		}
		return pj / l;
	}

	/**
	 * Computes the TPR using the top method.
	 * @param predResults The prediction results.
	 * @param results Iterator to write results.
	 */
	virtual void evalUsingTop(std::shared_ptr<PredResultsT> &predResults,
			PerfResults &results) {

		logger(logDebug, "Computing TopPrec using the method top...")

		auto remLinksMap = predResults->getTestData().getRemLinksMap();
		std::size_t fn = remLinksMap->size();
		if (l > fn) {
			std::cerr
					<< "# Warning in getTPR: l is larger than the number of positive links. The TPR cannot reach 1 in this case!\n";
		}

		auto net = predResults->getTestData().getObsNet();

		predResults->compTopScores(l);
		auto topEdgesBegin = predResults->topEdgesBegin();
		auto topEdgesEnd = predResults->topEdgesEnd();

//		std::cout << "#TPR top edges\n";
//		for (auto it = topEdgesBegin; it != topEdgesEnd; ++it) {
//			std::cout << net->start(*it) << "\t" << net->getLabel(net->end(*it))
//					<< std::endl;
//		}
//		std::cout << "-----------------\n";
//		std::cout << "#TPR removed edges\n";
//		for (auto it = remLinksMap->begin(); it != remLinksMap->end(); ++it) {
//			std::cout << net->getLabel(net->start(*it)) << "\t"
//					<< net->getLabel(net->end(*it)) << std::endl;
//		}
//		std::cout << "-----------------\n";
		std::size_t found = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for reduction (+:found) if (parallel)
#endif
		for (auto it = topEdgesBegin; it < topEdgesEnd; ++it) {
			if (remLinksMap->count(*it) > 0
					|| remLinksMap->count(net->reverseEdge(*it))) {
				found++;
			}
		}
		std::size_t k = predResults->getNbTop();
		std::size_t n = net->getNbNonEdges();
		results[name] = ((l - k) * (double) fn / n + found) / l;

		logger(logDebug, "Done")
	}

	/**
	 * Computes the TPR using the predict method.
	 * @param predResults The prediction results.
	 * @param results Iterator to write results.
	 */
	virtual void evalUsingPredict(std::shared_ptr<PredResultsT> &predResults,
			PerfResults &results) {

		logger(logDebug, "Computing TopPrec using the method predict...")

//		std::cout << "evalUsingPredict\n";
		predResults->compNegScores();
		predResults->sortNeg(SortOrder::Inc);
		auto negScoresBegin = predResults->negBegin();
		auto negScoresEnd = predResults->negEnd();

		predResults->compPosScores();
		predResults->sortPos(SortOrder::Inc);
		auto posScoresBegin = predResults->posBegin();
		auto posScoresEnd = predResults->posEnd();

//		Utilities::print(posScoresBegin, posScoresEnd, "Pos");
//		Utilities::print(negScoresBegin, negScoresEnd, "Neg");

		results[name] = getTPR(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd, l);

		logger(logDebug, "Done")
	}

	/**
	 * Computes the performance measure.
	 * @param predResults The prediction results.
	 * @param results Iterator to write results.
	 */
	virtual void eval(std::shared_ptr<PredResultsT> &predResults,
			PerfResults &results) {

		logger(logDebug, "Computing TopPrec...")

		if (useTopMethod) {
			evalUsingTop(predResults, results);
		} else {
			evalUsingPredict(predResults, results);
		}
		logger(logDebug, "Done")
	}

	/**
	 * Computes the area under the PR curve.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param posSortOrder The sorting order of positive scores. This may be modified by the method.
	 * @param negSortOrder The sorting order of negative scores. This may be modified by the method.
	 * @param results To write results.
	 */
	virtual void eval(ScoresIteratorT posScoresBegin,
			ScoresIteratorT posScoresEnd, ScoresIteratorT negScoresBegin,
			ScoresIteratorT negScoresEnd, SortOrder &posSortOrder,
			SortOrder &negSortOrder, PerfResults &results) {

		logger(logDebug, "Computing TopPrec...")
		if (posSortOrder != Inc) {
			Utilities::sort(posScoresBegin, posScoresEnd, SortOrder::Inc);
			posSortOrder = Inc;
		}
		if (negSortOrder != Inc) {
			Utilities::sort(negScoresBegin, negScoresEnd, SortOrder::Inc);
			negSortOrder = Inc;
		}

		results[name] = getTPR(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd, l);

		logger(logDebug, "Done")
	}

	/**
	 * @return True if the method top is used to compute TPR, false otherwise.
	 */
	bool isUseTopMethod() const {
		return useTopMethod;
	}

	/**
	 * Enable/disable the use of the method top is used to compute TPR.
	 * @param useTopMethod If true, enable the use of the method top, otherwise disable it.
	 */
	void setUseTopMethod(bool useTopMethod) {
		this->useTopMethod = useTopMethod;
	}

	/**
	 * @return Whether the performance measure requires the generation of negative set. The default value is true.
	 */
	virtual bool requiresNeg() const {
		return !useTopMethod;
	}

	/**
	 * @return Whether the performance measure requires network shuffling.
	 */
	virtual bool requiresShuffling() const {
		return useTopMethod;
	}

	/**
	 * Destructor.
	 */
	virtual ~TPR() = default;

};

/**
 * @brief General performance curve.
 * @tparam PredResultsT The prediction results type.
 */
template<typename PredResultsT = PredResults<>> class GCurve: public PerfCurve<
		PredResultsT> {

public:
	using ScoresIteratorT = typename PerfMeasure<PredResultsT>::ScoresIteratorT; /**< Scores iterator type. */

protected:
	using PerfMeasure<PredResultsT>::name; /**< The name of the performance measure. */
#ifdef WITH_OPENMP
	using PerfMeasure<PredResultsT>::parallel; /**< Enable/disable parallelism. */
#endif

protected:
	typename PerfLambda::PerfLambdaT const xLambda; /**< Lambda for computing x-coordinates. */
	typename PerfLambda::PerfLambdaT const yLambda; /**< Lambda for computing y-coordinates. */

public:
	/**
	 * Constructor.
	 * @param xLambda A lambda to compute the x-coordinates.
	 * @param yLambda A lambda to compute the y-coordinates.
	 */
	GCurve(typename PerfLambda::PerfLambdaT const &xLambda,
			typename PerfLambda::PerfLambdaT const &yLambda) :
			PerfCurve<PredResultsT>("GCurve"), xLambda(xLambda), yLambda(
					yLambda) {
	}

	/**
	 * Constructor.
	 * @param xLambda A lambda to compute the x-coordinates.
	 * @param yLambda A lambda to compute the y-coordinates.
	 * @param name The name of the performance measure.
	 */
	GCurve(typename PerfLambda::PerfLambdaT const &xLambda,
			typename PerfLambda::PerfLambdaT const &yLambda, std::string name) :
			PerfCurve<PredResultsT>(name), xLambda(xLambda), yLambda(yLambda) {

	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	GCurve(GCurve const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	GCurve& operator =(GCurve const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	GCurve(GCurve &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	GCurve& operator =(GCurve &&that) = default;

	/**
	 * Computes the area under the curve.
	 * @param predResults The prediction results.
	 * @param results Iterator to write results.
	 */
	virtual void eval(std::shared_ptr<PredResultsT> &predResults,
			PerfResults &results) {

		logger(logDebug, "Computing GCurve-AUC...")
#ifdef WITH_OPENMP
		results[name] = getGCurveAuc(getCurve(predResults), parallel);
#else
		results[name] = getGCurveAuc(getCurve(predResults));
#endif
		logger(logDebug, "Done")
	}

	/**
	 * Computes the area under the curve.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param posSortOrder The sorting order of positive scores. This may be modified by the method.
	 * @param negSortOrder The sorting order of negative scores. This may be modified by the method.
	 * @param results To write results.
	 */
	virtual void eval(ScoresIteratorT posScoresBegin,
			ScoresIteratorT posScoresEnd, ScoresIteratorT negScoresBegin,
			ScoresIteratorT negScoresEnd, SortOrder &posSortOrder,
			SortOrder &negSortOrder, PerfResults &results) {

		logger(logDebug, "Computing GCurve-AUC...")

#ifdef WITH_OPENMP
		results[name] = getGCurveAuc(
				getCurve(posScoresBegin, posScoresEnd, negScoresBegin,
						negScoresEnd, posSortOrder, negSortOrder), parallel);
#else
		results[name] = getGCurveAuc(
				getCurve(posScoresBegin, posScoresEnd, negScoresBegin,
						negScoresEnd, posSortOrder, negSortOrder));
#endif
		logger(logDebug, "Done")
	}

	/**
	 * Compute the threshold of the GCurve curve. Ranges must be sorted in increasing order.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 */
	template<typename ScoresIterator> static std::vector<double> getThresholds(
			ScoresIterator posScoresBegin, ScoresIterator posScoresEnd,
			ScoresIterator negScoresBegin, ScoresIterator negScoresEnd) {

		std::size_t P = posScoresEnd - posScoresBegin;
		std::size_t N = negScoresEnd - negScoresBegin;

		std::vector<double> allScores;
		allScores.resize(P + N);
		std::merge(posScoresBegin, posScoresEnd, negScoresBegin, negScoresEnd,
				allScores.begin());

		std::vector<double> unique;
		for (auto it = allScores.cbegin(); it < allScores.cend();) {
			unique.push_back(*it);
			auto rit = it + 1;
			while ((rit < posScoresEnd) && (*rit == *it)) {
				++rit;
			}
			it = rit;
		}
		unique.push_back(std::numeric_limits<double>::infinity());
//		Utils::print(unique.begin(), unique.end(), "thresholds");
		return unique;
	}

	/**
	 * Compute the curve. Both ranges must be sorted.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param parallel Whether to run in parallel.
	 */
	template<typename ScoresIterator> std::vector<std::pair<double, double>> getGCurve(
			ScoresIterator posScoresBegin, ScoresIterator posScoresEnd,
			ScoresIterator negScoresBegin, ScoresIterator negScoresEnd,
			bool parallel = false) {

		std::size_t P = posScoresEnd - posScoresBegin;
		std::size_t N = negScoresEnd - negScoresBegin;

		auto threshs = getThresholds(posScoresBegin, posScoresEnd,
				negScoresBegin, negScoresEnd);

		std::vector<std::pair<double, double>> curve;
		curve.resize(threshs.size());
		std::size_t maxInd = threshs.size() - 1;

#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
		for (std::size_t i = 0; i < threshs.size(); i++) {

			auto itpl = std::lower_bound(posScoresBegin, posScoresEnd,
					threshs[i]);
			std::size_t tp = posScoresEnd - itpl;
			std::size_t fn = P - tp;

			auto itnl = std::lower_bound(negScoresBegin, negScoresEnd,
					threshs[i]);
			std::size_t fp = negScoresEnd - itnl;
			std::size_t tn = N - fp;
			curve[maxInd - i].first = xLambda(tp, fn, tn, fp, P, N);
			curve[maxInd - i].second = yLambda(tp, fn, tn, fp, P, N);
		}

		// Sort the coordinates
		std::sort(curve.begin(), curve.end());
		return curve;
	}

	/**
	 * Computes the performance curve.
	 * @param predResults The prediction results.
	 * @return A curve in the form of an std::vector of pairs representing the x and y coordinates.
	 */
	virtual std::vector<std::pair<double, double>> getCurve(
			std::shared_ptr<PredResultsT> &predResults) {

		predResults->compNegScores();
		predResults->sortNeg(SortOrder::Inc);
		auto negScoresBegin = predResults->negBegin();
		auto negScoresEnd = predResults->negEnd();

		predResults->compPosScores();
		predResults->sortPos(SortOrder::Inc);
		auto posScoresBegin = predResults->posBegin();
		auto posScoresEnd = predResults->posEnd();

#ifdef WITH_OPENMP
		return getGCurve(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd, parallel);
#else
		return getGCurve(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd);
#endif
	}

	/**
	 * Computes the performance curve.
	 * @param posScoresBegin Iterator to the first positive score.
	 * @param posScoresEnd Iterator to one-past-the-last positive score.
	 * @param negScoresBegin Iterator to the first negative score.
	 * @param negScoresEnd Iterator to one-past-the-last negative score.
	 * @param posSortOrder The sorting order of positive scores. This may be modified by the method.
	 * @param negSortOrder The sorting order of negative scores. This may be modified by the method.
	 * @return A curve in the form of an std::vector of pairs representing the x and y coordinates.
	 */
	virtual std::vector<std::pair<double, double>> getCurve(
			ScoresIteratorT posScoresBegin, ScoresIteratorT posScoresEnd,
			ScoresIteratorT negScoresBegin, ScoresIteratorT negScoresEnd,
			SortOrder &posSortOrder, SortOrder &negSortOrder) {

		if (posSortOrder != Inc) {
			Utilities::sort(posScoresBegin, posScoresEnd, SortOrder::Inc);
			posSortOrder = Inc;
		}
		if (negSortOrder != Inc) {
			Utilities::sort(negScoresBegin, negScoresEnd, SortOrder::Inc);
			negSortOrder = Inc;
		}

#ifdef WITH_OPENMP
		return getGCurve(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd, parallel);
#else
		return getGCurve(posScoresBegin, posScoresEnd, negScoresBegin,
				negScoresEnd);
#endif
	}

	/**
	 * Computes the area under the curve using integration by linear interpolation.
	 * @param curve The curve represented by its x-y coordinates.
	 * @param parallel Whether to run in parallel.
	 * @return The area under the curve.
	 */
	static double getGCurveAuc(
			std::vector<std::pair<double, double>> const &curve, bool parallel =
					false) {

		logger(logDebug, "Computing AUC")

		double auc = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for reduction(+ : auc) if (parallel)
#endif
		for (std::size_t i = 0; i < curve.size() - 1; i++) {
			auc += 0.5 * (curve[i + 1].first - curve[i].first)
					* (curve[i].second + curve[i + 1].second);
		}

		logger(logDebug, "Done computing AUC")
		return auc;
	}

	/**
	 * Destructor.
	 */
	virtual ~GCurve() = default;
};

} /* namespace LinkPred */

#endif /* INCLUDE_PERFMEASURE_HPP_ */
