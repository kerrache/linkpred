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
 * @brief Includes the headers related to core classes.
 */

#ifndef PERFEVALUATOR_HPP_
#define PERFEVALUATOR_HPP_

#include <linkpred/predictors/ulpredictor.hpp>
#include "linkpred/predictors/dlpredictor.hpp"
#include "linkpred/perf/predresults.hpp"
#include "linkpred/perf/perfmeasure.hpp"
#include "linkpred/perf/networkmanipulator.hpp"
#include "linkpred/utils/randomgen.hpp"
#include "linkpred/utils/log.hpp"
#include "linkpred/utils/utilities.hpp"
#include <string>
#include <memory>
#include <vector>
#include <chrono>
#include <map>
#include <stdexcept>
#include <iostream>
#include <fstream>

namespace LinkPred {

/**
 * @brief Performance evaluator.
 * @tparam TestDataT The test data type.
 * @tparam PredResultsT The prediction results type.
 * @tparam PerfMeasureT The performance measure type.
 */
template<typename TestDataT = TestData<>, typename LPredictorT = ULPredictor<>,
		typename PredResultsT = PredResults<TestDataT, LPredictorT>,
		typename PerfMeasureT = PerfMeasure<PredResultsT>> class PerfEvaluator {

protected:
	std::vector<std::shared_ptr<LPredictorT>> predictors; /**< To store predictors. */
	std::vector<std::shared_ptr<PerfMeasureT>> measures; /**< To store performance measures. */
	TestDataT testData; /**< The test data. */
#ifdef WITH_OPENMP
	bool parallel = false; /**< To enable/disable parallelism. */
#endif
	bool timingEnabled = false; /**< To enable/disable timing. */
	std::map<std::string, double> results; /**< Performance results. */

public:
	/**
	 * Constructor.
	 * @param testData The test data.
	 */
	PerfEvaluator(TestDataT testData) :
			testData(testData) {
		if (!testData.isLocked()) {
			throw std::runtime_error("TestData must be locked");
		}
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	PerfEvaluator(PerfEvaluator const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	PerfEvaluator & operator =(PerfEvaluator const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	PerfEvaluator(PerfEvaluator && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	PerfEvaluator & operator =(PerfEvaluator && that) = default;

	/**
	 * Add a predictor.
	 * @param predictor A link predictor.
	 * @return An ID.
	 */
	std::size_t addPredictor(std::shared_ptr<LPredictorT> predictor) {
		predictors.push_back(predictor);
		return predictors.size() - 1;
	}

	/**
	 * Add a performance measure.
	 * @param measure A performanmce measure.
	 * @return An ID.
	 */
	std::size_t addPerfMeasure(std::shared_ptr<PerfMeasureT> measure) {
		measures.push_back(measure);
		return measures.size() - 1;
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

	/**
	 * @return Ther number of predictors.
	 */
	std::size_t getNbPredictors() const {
		return predictors.size();
	}

	/**
	 * @return The number of performance measures.
	 */
	std::size_t getNbPerfMeasures() const {
		return measures.size();
	}

	/**
	 * @param i Index of the predictor.
	 * @return The predictor of index i.
	 */
	auto getPredictor(std::size_t i) const {
		return predictors[i];
	}

	/**
	 * Run the performance evaluation.
	 */
	void eval() {

		logger(logDebug, "Computing performance measures...")
		if (timingEnabled) {
			evalTiming();
		} else {
			evalNoTiming();
		}
		logger(logDebug, "Done")
	}

	/**
	 * Runs the performance evaluation without timing.
	 */
	void evalNoTiming();

	/**
	 * Runs the performance evaluation with timing.
	 */
	void evalTiming();

	/**
	 * @return Whether timing is enabled.
	 */
	bool isTimingEnabled() const {
		return timingEnabled;
	}

	/**
	 * Enable/disable timing.
	 * @param timingEnabled The new value.
	 */
	void setTimingEnabled(bool timingEnabled) {
		this->timingEnabled = timingEnabled;
	}

	/**
	 * @return An iterator to the first performance result.
	 */
	auto resultsBegin() const {
		return results.cbegin();
	}

	/**
	 * @return An iterator to one-past-the-last first performance result.
	 */
	auto resultsEnd() const {
		return results.cend();
	}

	/**
	 * Destructor.
	 */
	virtual ~PerfEvaluator() = default;
};

/**
 * @brief Structure storing experiment description.
 */
template<typename NetworkT = UNetwork<>> struct PerfeEvalExpDescp {

#ifdef WITH_OPENMP
	bool parTestRuns =  false; /**< To enable/disable shared-memory parallelism over test runs. */
	bool parPredictors = false; /**< Whether to run different predictors in parallel. */
#endif
#ifdef WITH_MPI
	bool distributed = false; /**< Whether the experiment is run distributively. */
	MPI_Comm comm = MPI_COMM_WORLD; /**< The MPI communicator. */
#endif
	bool timingEnabled = false; /**< Enable/disable timing. */
	std::ostream* out = &std::cout; /**< Output file. */
	std::size_t nbTestRuns = 1; /**< Number of test runs. */
	std::shared_ptr<NetworkT> refNet; /**< Reference network. */
	bool keepConnected = false; /**< Whether to keep the network connected. */
	double fnRatio = 1; /**< Ratio of false negatives used in the test set. */
	double tnRatio = 1; /**< Ratio of true negatives used in the test set. */
	double ratioStart = 0.1; /**< Start value of the ratio of removed edges. */
	double ratioEnd = 0.1; /**< End value of the ratio of removed edges. This is adjusted if keepConnected is set to true. */
	double ratioStep = 0.1; /**< Step size of the ratio of removed edges. */
	long int seed = 0; /**< Seed for the random number generator. */
};

/**
 * @brief Factory class to create link predictors and performance measures.
 * @tparam NetworkT The network data type.
 * @tparam LPredictorT The link predictor type.
 * @tparam TestDataT The test data type.
 * @tparam PerfMeasureT The performance measure type.
 */
template<typename NetworkT = UNetwork<>, typename LPredictorT = ULPredictor<>,
		typename TestDataT = TestData<>, typename PerfMeasureT = PerfMeasure<
				PredResults<TestDataT, LPredictorT>> > class PEFactory {
public:
	/**
	 * @return A vector of predictors.
	 */
	virtual std::vector<std::shared_ptr<LPredictorT>> getPredictors(
			std::shared_ptr<NetworkT const> obsNet) = 0;

	/**
	 * @return A vector of performance measures.
	 */
	virtual std::vector<std::shared_ptr<PerfMeasureT>> getPerfMeasures(
			TestDataT const & testData) =0;

	/**
	 * Destructor.
	 */
	virtual ~PEFactory() = default;
};

/**
 * @brief Performance evaluation experiment.
 * @tparam NetworkT The network data type.
 * @tparam TestDataT The test data type.
 * @tparam PredResultsT The prediction results type.
 * @tparam PerfMeasureT The performance measure type.
 */
template<typename NetworkT = UNetwork<>, typename TestDataT = TestData<>,
		typename LPredictorT = ULPredictor<>,
		typename PredResultsT = PredResults<TestDataT, LPredictorT>,
		typename PerfMeasureT = PerfMeasure<PredResultsT>> class PerfEvalExp {

protected:
	PerfeEvalExpDescp<NetworkT> ped; /**< Description of the experiment. */
	std::shared_ptr<PEFactory<NetworkT, LPredictorT, TestDataT, PerfMeasureT>> factory; /**< Factory to create link predictors and performance measures. */
	std::vector<std::map<std::string, double>> results; /**< Performance results. */
	int outPrec = 4; /**< Precision of print (for double). */

public:
	/**
	 * Constructor.
	 * @param testData The test data.
	 */
	PerfEvalExp(PerfeEvalExpDescp<NetworkT> const & ped,
			std::shared_ptr<
					PEFactory<NetworkT, LPredictorT, TestDataT, PerfMeasureT>> const & factory) :
			ped(ped), factory(factory) {
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	PerfEvalExp(PerfEvalExp const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	PerfEvalExp & operator =(PerfEvalExp const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	PerfEvalExp(PerfEvalExp && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	PerfEvalExp & operator =(PerfEvalExp && that) = default;

	/**
	 * Run the performance evaluation experiment.
	 */
	void run();

	/**
	 * Runs the performance evaluation without timing.
	 */
	void runNoTiming();

	/**
	 * Runs the performance evaluation with timing.
	 */
	void runTiming();

	/**
	 * @return An iterator to the first performance result.
	 */
	auto resultsBegin() const {
		return results.cbegin();
	}

	/**
	 * @return An iterator to one-past-the-last first performance result.
	 */
	auto resultsEnd() const {
		return results.cend();
	}

	/**
	 * @return The fatory.
	 */
	const std::shared_ptr<
			PEFactory<NetworkT, LPredictorT, TestDataT, PerfMeasureT> >& getFactory() const {
		return factory;
	}

	/**
	 * @return Output precision (for double).
	 */
	int getOutPrec() const {
		return outPrec;
	}

	/**
	 * @param outPrec Output precision (for double).
	 */
	void setOutPrec(int outPrec) {
		this->outPrec = outPrec;
	}

	/**
	 * @return The experience descriptor.
	 */
	const PerfeEvalExpDescp<NetworkT>& getPed() const {
		return ped;
	}

	/**
	 * Destructor.
	 */
	virtual ~PerfEvalExp() = default;

};

}
/* namespace LinkPred */

#endif /* PERFEVALUATOR_HPP_ */
