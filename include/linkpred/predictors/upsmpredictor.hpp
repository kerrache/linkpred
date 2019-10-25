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
 * @brief Contains the implementation of a popularity similarity link predictor using MDS.
 */

#ifndef UPSMPREDICTOR_HPP_
#define UPSMPREDICTOR_HPP_
#include <linkpred/predictors/ulpredictor.hpp>
#include "linkpred/utils/randomgen.hpp"
#include <cmath>

namespace LinkPred {
/**
 * @brief Popularity-similarity link predictor using MDS.
 * @tparam NetworkT The network type.
 * @tparam EdgesRandomIteratorT A random iterator type used to iterate on edges.
 * @tparam ScoresRandomIteratorT A random iterator type used to iterate on scores.
 */
template<typename NetworkT = UNetwork<>,
		typename EdgesRandomIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::const_iterator,
		typename ScoresRandomIteratorT = typename std::vector<double>::iterator,
		typename EdgesRandomOutputIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::iterator> class UPSMPredictor: public ULPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT> {

	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::net; /**< The network. */
#ifdef WITH_OPENMP
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::parallel; /**< Whether the predictor runs in parallel. */
#endif
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::name; /**< The name of the predictor. */
	using NodeIdType = typename ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::NodeIdType; /**< The node IDs type. */
	using EdgeType = typename ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::EdgeType; /**< The edges type. */

	/**
	 * @brief Enumeration of the possible link statuses
	 */
	enum LinkStatus {
		Pos, /**< Positive link. */
		Neg, /**< Negative link. */
		Unk, /**< Unknown. */
	};

protected:
	long int seed; /**< Random generator seed. */
	RandomGen rng; /**< Random number generator. */
	std::size_t dim = 0; /**< The dimension. */
	std::vector<double> coords; /**< The coordinates of the nodes. */
	std::vector<double> logDeg; /**< Log of the observed degrees of the nodes. */
	std::vector<double> logCDeg; /**< Log of the computed degrees of the nodes. */
	std::vector<double> logCDegProd; /**< Log of the product of the computed degrees. */
	std::vector<double> sqDist; /**< Squared distances between nodes. */
	std::vector<LinkStatus> linkStatus; /**< Link statuses. */
	std::vector<double> weight; /**< Link weight in MDS. */
	std::vector<double> probs; /**< Probabilities of connection. */
	std::size_t nbNodes = 0; /**< Number of nodes. */
	std::size_t nbKnownCouples = 0; /**< Number of known couples. */
	std::size_t nbCouples = 0; /**< Number of couples. */
	double wPos = -1; /**< Weight of positive links. */
	double wNeg = -1; /**< Weight of negative links. */
	std::size_t nbMdsRuns = 1; /**< Number of times MDS is run. */
	double globObj = 0; /**< Global objective. */
	double oldGlobObj = std::numeric_limits<double>::max(); /**< Old objective value. */
	double deltaGlobObjTol = 1.0e-7; /**< Tolerance in detltaGlobObj. */
	std::size_t maxIter = 100; /**< Maximum number of iterations. */
	std::size_t iter = 0; /**< The number of iterations. */
	double cl = 0; /**< 1 - confidence level. */
	double sumW = 0; /**< Sum of weights. */
	double eta = 3; /**< Factor between connected and disconnected links. */
	std::size_t maxDim = 5; /**<Maximum dimension to try. */

	/**
	 * Termination test.
	 * @return True if program solver must terminate.
	 */
	bool terminate() const;

	/**
	 * Compute global objective: Error between observed and computed probabilities.
	 */
	void compGlobObj();

	/**
	 * Detect dimension.
	 */
	void detectDim();

	/**
	 * Compute the distances.
	 */
	void compDist();

	/**
	 * Solve the MDS problem.
	 * @param init If true the starting point is initialized by the method.
	 * @return The MDS objective.
	 */
	double mds(bool init);

	/**
	 * Estimate the degrees.
	 * @param init If true the starting point is initialized by the method.
	 */
	void compDeg(bool init);

	/**
	 * Remove outliers.
	 */
	void removeOutliers();

	/**
	 * The PDF of the sum of two IID exponentials.
	 * @param x The point at which the PDF is computed.
	 * @param lm The parameter of the exponential.
	 */
	inline double sumExpPdf(double x, double lm) {
		return lm * lm * x * std::exp(-lm * x);
	}

public:

	/**
	 * @param net The network.
	 * @param seed The random number generator's seed.
	 */
	UPSMPredictor(std::shared_ptr<NetworkT const> net, long int seed) :
			ULPredictor<NetworkT, EdgesRandomIteratorT>(net), seed(seed), rng(
					seed) {
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	UPSMPredictor(UPSMPredictor const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	UPSMPredictor & operator =(UPSMPredictor const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	UPSMPredictor(UPSMPredictor && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	UPSMPredictor & operator =(UPSMPredictor && that) = default;

	/**
	 * Initialize the solver.
	 */
	virtual void init();

	/**
	 * Learning.
	 */
	virtual void learn();

	/**
	 * Compute the score of a single edge.
	 * @param e The edge.
	 * @return The score of e.
	 */
	virtual double score(EdgeType const & e) {
		auto k = net->coupleOrd(e);
		return 1.0 / (1 + std::sqrt(sqDist[k]) / std::exp(logCDegProd[k]));
	}

	/**
	 * Predict the links.
	 * @param begin Beginning of the links to be predicted.
	 * @param end end of the links to be predicted.
	 * @param scores Beginning of scores.
	 */
	virtual void predict(EdgesRandomIteratorT begin, EdgesRandomIteratorT end,
			ScoresRandomIteratorT scores);

	/**
	 * @return The dimension of the embedding space.
	 */
	std::size_t getDim() const {
		return dim;
	}

	/**
	 * Set the dimension of the embedding space.
	 * @param dim The new dimension of the embedding space.
	 */
	void setDim(std::size_t dim = 2) {
		this->dim = dim;
	}

	/**
	 * @return The weight of negative links.
	 */
	double getWNeg() const {
		return wNeg;
	}

	/**
	 * Set the weight of negative links.
	 * @param wNeg The new weight of negative links.
	 */
	void setWNeg(double wNeg) {
		this->wNeg = wNeg;
	}

	/**
	 * @return The weight of positive links.
	 */
	double getWPos() const {
		return wPos;
	}

	/**
	 * Set the weight of positive links.
	 * @param wPos The new weight of positive links.
	 */
	void setWPos(double wPos) {
		this->wPos = wPos;
	}

	/**
	 * @return The number of runs used in the initial multidimensional scaling.
	 */
	std::size_t getNbMdsRuns() const {
		return nbMdsRuns;
	}

	/**
	 * Set the number of runs used in the initial multidimensional scaling.
	 * @param nbMdsRuns The new number of runs used in the initial multidimensional scaling.
	 */
	void setNbMdsRuns(std::size_t nbMdsRuns) {
		this->nbMdsRuns = nbMdsRuns;
	}

	/**
	 * @return The tolerance in difference in the global objective value.
	 */
	double getDeltaGlobObjTol() const {
		return deltaGlobObjTol;
	}

	/**
	 * Set the tolerance in difference in the global objective value.
	 * @param deltaGlobObjTol The new tolerance in difference in the global objective value.
	 */
	void setDeltaGlobObjTol(double deltaGlobObjTol) {
		this->deltaGlobObjTol = deltaGlobObjTol;
	}

	/**
	 * @return The number of iterations.
	 */
	std::size_t getIter() const {
		return iter;
	}

	/**
	 * @return The maximum number of iterations.
	 */
	std::size_t getMaxIter() const {
		return maxIter;
	}

	/**
	 * Set the maximum number of iterations.
	 * @param maxIter The new maximum number of iterations.
	 */
	void setMaxIter(std::size_t maxIter) {
		this->maxIter = maxIter;
	}

	/**
	 * @return The global objective value.
	 */
	double getGlobObj() const {
		return globObj;
	}

	/**
	 * @return The confidence level.
	 */
	double getCl() const {
		return cl;
	}

	/**
	 * Set the confidence level.
	 * @param cl The new confidence level.
	 */
	void setCl(double cl) {
		this->cl = cl;
	}

	/**
	 * @return The parameter eta.
	 */
	double getEta() const {
		return eta;
	}

	/**
	 * Set the parameter eta.
	 * @param eta The new value of the parameter eta.
	 */
	void setEta(double eta) {
		this->eta = eta;
	}

	/**
	 * @return The maximum dimension.
	 */
	std::size_t getMaxDim() const {
		return maxDim;
	}

	/**
	 * Set the maximum dimension.
	 * @param maxDim The new maximum dimension.
	 */
	void setMaxDim(std::size_t maxDim) {
		this->maxDim = maxDim;
	}

	/**
	 * Destructor.
	 */
	virtual ~UPSMPredictor() = default;
};

}
/* namespace LinkPred */

#endif /* UPSMPREDICTOR_HPP_ */
