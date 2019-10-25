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
 * @brief Contains the implementation of a scalable popularity similarity link predictor.
 */

#ifndef UMPSPREDICTOR_HPP_
#define UMPSPREDICTOR_HPP_

#include <linkpred/predictors/ulpredictor.hpp>
#include "linkpred/core/netdistcalculator.hpp"
#include "linkpred/perf/perfmeasure.hpp"
#include "linkpred/utils/log.hpp"
#include <memory>
#include <cmath>
#include <limits>

namespace LinkPred {

/**
 * @brief A scalable popularity similarity link predictor.
 * @tparam NetworkT The network type.
 * @tparam EdgesRandomIteratorT A random iterator type used to iterate on edges.
 * @tparam ScoresRandomIteratorT A random iterator type used to iterate on scores.
 */
template<typename NetworkT = UNetwork<>,
		typename EdgesRandomIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::const_iterator,
		typename ScoresRandomIteratorT = typename std::vector<double>::iterator,
		typename EdgesRandomOutputIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::iterator> class UMPSPredictor: public ULPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT> {

	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::net;
#ifdef WITH_OPENMP
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::parallel;
#endif
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::name;
	using NodeIdType = typename ULPredictor<NetworkT, EdgesRandomIteratorT,ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::NodeIdType;
	using EdgeType = typename ULPredictor<NetworkT, EdgesRandomIteratorT,ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::EdgeType;

public:
	/**
	 * @brief An enumeration of the different landmark positioning strategies.
	 */
	enum LandmarkStrategy {
		Random, /**< Landmarks are chosen randomly. */
		Hub, /**< The nodes with the highest degree are chosen. */
		IHub /**< The nodes with the lowest degree are chosen. */
	};

	/**
	 * @brief An enumeration of different methods to find lambda.
	 */
	enum LambdaMethod {
		User, /**< Use the value fixed by the user. */
		MeanApp, /**< Lambda is approximated by an average. */
		Scan, /**< Scan specified values and choose the best. */
		Opt /**< Find lambda by optimization. */
	};

	/**
	 * @brief An enumeration of edge length methods.
	 */
	enum EdgeLengthMethod {
		PAT, /**< Degree products. */
		RA1, /**< RA1 method. */
		RA2 /**< RA2 method. */
	};

protected:
	std::shared_ptr<NetworkT const> oriNet; /**< Original network. */
	std::shared_ptr<NetworkT> modNet; /**< Network obtained after removing links. */
	double remRatio; /**< Ratio of removed links. */
	double lambda = 0.5; /**< The parameter lambda. */
	LambdaMethod lambdaMethod = User; /**< How lambda is estimated. */
	double lambdaPosEstRatio = 0.1; /**< Ratio of positive links used to estimate lambda (used only if lambdaMethod != User).*/
	double lambdaNegEstRatio = 1; /**< Ratio of negative links used to estimate lambda (used only if lambdaMethod != User).*/
	double negScore = 0.01; /**< Score given to a negative edge. Must be in ]0, 1[. */
	double posScore = 0.99; /**< Score given to a positive edge. Must be in ]0, 1[. */
	std::shared_ptr<PerfMeasure<>> perfMeasure; /**< The performance measure used to estimate lambda (used only if estimateLambda is true).*/
	double minLambda = 0.1; /**< Minimum lambda value to try. Used when betMethod is set to Scan. */
	double maxLambda = 0.9; /**< Maximum lambda value to try. Used when betMethod is set to Scan. */
	double lambdaStep = 0.1; /**< The step size when scanning for lambda. Used when betMethod is set to Scan. */
	double tol = 1e-5; /**< Tolerance when optimizing for bet. Used only when lambdaMethod == Opt. */
	long int seed = 0; /**< The random number generator seed. */
	RandomGen rng; /**< The random number generator. */
	Dijkstra<NetworkT, double, std::size_t> dijkstra; /**< Dijkstra's algorithm. */
	typename NetDistCalculator<NetworkT, double, std::size_t>::EdgeLengthMapSP length; /**< The length map. */
	std::shared_ptr<NetDistCalculator<NetworkT, double, std::size_t>> distCalc; /**< The distance calculator. */
	CacheLevel cacheLevel = NodeCache; /**< The cache level. */
	bool asp = false; /**< Whether to use approximate shortest paths instead of exact ones. */
	double landmarkRatio = 0; /**< The ratio of nodes used as landmarks if asp is on. */
	LandmarkStrategy landmarkStrategy = Random; /**< Landmark positioning strategy if asp is on. */
	std::set<NodeIdType> landmarks; /**< The landmarks. */
	EdgeLengthMethod edgeLengthMethod = PAT; /**< The method used to compute edge lengths. */

	/**
	 * Compute the length of the edge according to the PAT method.
	 * @param edge The edge.
	 * @return The edge length using PAT method.
	 */
	inline double PATEdgeLength(EdgeType const & edge) const {
		double ki = net->getDeg(NetworkT::start(edge));
		double kj = net->getDeg(NetworkT::end(edge));
		return ki * kj + std::numeric_limits<double>::epsilon();
	}

	/**
	 * Compute the length of the edge according to the RA1 method.
	 * @param edge The edge.
	 * @return The edge length using RA1 method.
	 */
	inline double RA1EdgeLength(EdgeType const & edge) const {
		auto srcNode = NetworkT::start(edge);
		auto endNode = NetworkT::end(edge);
		double ki = net->getDeg(srcNode);
		double kj = net->getDeg(endNode);
		double cn = net->getNbCommonNeighbors(srcNode, endNode);
		return (ki + kj + ki * kj) / (1 + cn)
				+ std::numeric_limits<double>::epsilon();
	}

	/**
	 * Compute the length of the edge according to the RA1 method.
	 * @param edge The edge.
	 * @return The edge length using RA1 method.
	 */
	inline double RA2EdgeLength(EdgeType const & edge) const {
		auto srcNode = NetworkT::start(edge);
		auto endNode = NetworkT::end(edge);
		double cn = net->getNbCommonNeighbors(srcNode, endNode);
		double eij;
		if (net->isEdge(srcNode, endNode)) {
			eij = 1;
		} else {
			eij = 0;
		}
		double ei = net->getDeg(srcNode) - cn - eij;
		double ej = net->getDeg(endNode) - cn - eij;
		return (1 + ei + ej + ei * ej) / (1 + cn)
				+ std::numeric_limits<double>::epsilon();
	}

	/**
	 * Create the length map from network topology.
	 */
	void createLengthMap();

	/**
	 * Position the landmarks.
	 */
	void setLandmarks();

	/**
	 * Estimate lambda from data.
	 */
	void estimateLambda();

	/**
	 * Find lambda by approximation.
	 */
	void lambdaMeanApp();

	/**
	 * Find lambda by scanning.
	 */
	void lambdaScan();

	/**
	 * Find lambda by optimization.
	 */
	void lambdaOpt();

public:
	/**
	 * @param net The network.
	 * @param seed The random number generator's seed.
	 */
	UMPSPredictor(std::shared_ptr<NetworkT const> net, double remRatio,
			long int seed) :
			ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
					EdgesRandomOutputIteratorT>(net), oriNet(net), remRatio(
					remRatio), seed(seed), rng(seed) {
		name = "MPS";

		// Copy of the network.
		modNet.reset(new NetworkT());
		for (auto it = net->nodesBegin(); it != net->nodesEnd(); ++it) {
			modNet->addNode(it->second);
		}

		for (auto it = net->rndEdgesBegin(1 - remRatio, rng.getInt());
				it != net->rndEdgesEnd(); ++it) {
			modNet->addEdge(net->start(*it), net->end(*it));
		}
		modNet->assemble();
		net = modNet;
		dijkstra.setNet(net);
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	UMPSPredictor(UMPSPredictor const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	UMPSPredictor & operator =(UMPSPredictor const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	UMPSPredictor(UMPSPredictor && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	UMPSPredictor & operator =(UMPSPredictor && that) = default;

	/**
	 * Initialize the predictor.
	 */
	virtual void init();

	/**
	 * Learn.
	 */
	virtual void learn();

	/**
	 * Compute the score of a single edge.
	 * @param e The edge.
	 * @return The score of e.
	 */
	virtual double score(EdgeType const & e) {

		double sc = 0;
		auto srcNode = NetworkT::start(e);
		auto endNode = NetworkT::end(e);
		switch (edgeLengthMethod) {
		case PAT: {
			if (net->isEdge(srcNode, endNode)) {
				auto res = distCalc->getIndDist(srcNode, endNode);
				sc =
						1.0
								/ (1
										+ std::pow(res.first, lambda)
												* std::pow(PATEdgeLength(e),
														lambda - 1));
			} else {
				auto res = distCalc->getDist(srcNode, endNode);
				sc =
						1.0
								/ (1
										+ std::pow(res.first, lambda)
												* std::pow(PATEdgeLength(e),
														lambda - 1));
			}
			break;
		}

		case RA1: {
			if (net->isEdge(srcNode, endNode)) {
				auto res = distCalc->getIndDist(srcNode, endNode);
				sc =
						1.0
								/ (1
										+ std::pow(res.first, lambda)
												* std::pow(RA1EdgeLength(e),
														lambda - 1));
			} else {
				auto res = distCalc->getDist(srcNode, endNode);
				sc =
						1.0
								/ (1
										+ std::pow(res.first, lambda)
												* std::pow(RA1EdgeLength(e),
														lambda - 1));
			}

			break;
		}
		case RA2: {
			if (net->isEdge(srcNode, endNode)) {
				auto res = distCalc->getIndDist(srcNode, endNode);
				sc =
						1.0
								/ (1
										+ std::pow(res.first, lambda)
												* std::pow(RA2EdgeLength(e),
														lambda - 1));
			} else {
				auto res = distCalc->getDist(srcNode, endNode);
				sc =
						1.0
								/ (1
										+ std::pow(res.first, lambda)
												* std::pow(RA2EdgeLength(e),
														lambda - 1));
			}

			break;
		}

		default:
			throw std::invalid_argument("Unknown edge length method");
		}

		return sc;
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
	 * @return The length map.
	 */
	auto getLength() const {
		return length;
	}

	/**
	 * @return Whether approximate shortest path distances are used.
	 */
	bool getAsp() const {
		return asp;
	}

	/**
	 * Set Whether approximate shortest path distances are used.
	 * @param asp Whether approximate shortest path distances should be used.
	 */
	void setAsp(bool asp) {
		this->asp = asp;
	}

	/**
	 * @return The landmark ratio.
	 */
	double getLandmarkRatio() const {
		return landmarkRatio;
	}

	/**
	 * Set the landmark ratio.
	 * @param landmarkRatio The new landmark ratio.
	 */
	void setLandmarkRatio(double landmarkRatio) {
		if (landmarkRatio < 0 || landmarkRatio > 1) {
			throw std::invalid_argument(
					"Invalid landmark ratio, must be between 0 and 1 inclusive");
		}
		this->landmarkRatio = landmarkRatio;
	}

	/**
	 * @return The landmark strategy.
	 */
	LandmarkStrategy getLandmarkStrategy() const {
		return landmarkStrategy;
	}

	/**
	 * Set the landmark strategy.
	 * @param landmarkStrategy The new landmark strategy.
	 */
	void setLandmarkStrategy(LandmarkStrategy landmarkStrategy) {
		this->landmarkStrategy = landmarkStrategy;
	}

	/**
	 * @return The parameter lambda.
	 */
	double getLambda() const {
		return lambda;
	}

	/**
	 * Set the parameter lambda.
	 * @param lambda The new value of the parameter lambda.
	 */
	void setLambda(double lambda) {
		this->lambda = lambda;
	}

	/**
	 * @return The distances cache level.
	 */
	CacheLevel getCacheLevel() const {
		return cacheLevel;
	}

	/**
	 * Set the distances cache level.
	 * @param cacheLevel  The new distances cache level.
	 */
	void setCacheLevel(CacheLevel cacheLevel) {
		this->cacheLevel = cacheLevel;
	}

	/**
	 * @return Score given to a negative edge.
	 */
	double getNegScore() const {
		return negScore;
	}

	/**
	 * Set the score given to a negative edge.
	 * @param negScore The new score given to a negative edge.
	 */
	void setNegScore(double negScore) {
		if (negScore >= 1 || negScore <= 0) {
			throw std::out_of_range("negScore must be in ]0, 1[.");
		}
		this->negScore = negScore;
	}

	/**
	 * @return Score given to a positive edge.
	 */
	double getPosScore() const {
		return posScore;
	}

	/**
	 * Set the score given to a positive edge.
	 * @param posScore The new score given to a positive edge.
	 */
	void setPosScore(double posScore) {
		if (posScore >= 1 || posScore <= 0) {
			throw std::out_of_range("posScore must be in ]0, 1[.");
		}
		this->posScore = posScore;
	}

	/**
	 * @return The lambda estimation method.
	 */
	LambdaMethod getLambdaMethod() const {
		return lambdaMethod;
	}

	/**
	 * Set the lambda estimation method.
	 * @param lambdaMethod The new lambda estimation method.
	 */
	void setLambdaMethod(LambdaMethod lambdaMethod) {
		this->lambdaMethod = lambdaMethod;
	}

	/**
	 * @return The step size when scanning for lambda. Used when betMethod is set to Scan.
	 */
	double getLambdaStep() const {
		return lambdaStep;
	}

	/**
	 * Set the step size when scanning for lambda. Used when betMethod is set to Scan.
	 * @param lambdaStep The new step size when scanning for lambda. Used when betMethod is set to Scan.
	 */
	void setLambdaStep(double lambdaStep) {
		this->lambdaStep = lambdaStep;
	}

	/**
	 * @return The maximum lambda value to try. Used when betMethod is set to Scan.
	 */
	double getMaxLambda() const {
		return maxLambda;
	}

	/**
	 * Set the maximum lambda value to try. Used when betMethod is set to Scan.
	 * @param maxLambda The new maximum lambda value to try. Used when betMethod is set to Scan.
	 */
	void setMaxLambda(double maxLambda) {
		this->maxLambda = maxLambda;
	}

	/**
	 * @return The minimum lambda value to try. Used when betMethod is set to Scan.
	 */
	double getMinLambda() const {
		return minLambda;
	}

	/**
	 * Set the minimum lambda value to try. Used when betMethod is set to Scan.
	 * @param minLambda The new minimum lambda value to try. Used when betMethod is set to Scan.
	 */
	void setMinLambda(double minLambda) {
		this->minLambda = minLambda;
	}

	/**
	 * @return The performance measure used to determine the best edge score method.
	 */
	auto getPerfMeasure() const {
		return perfMeasure;
	}

	/**
	 * Set the performance measure used to determine the best edge score method.
	 * @param perfMeasure The new performance measure used to determine the best edge score method.
	 */
	void setPerfMeasure(std::shared_ptr<PerfMeasure<>> perfMeasure) {
		this->perfMeasure = perfMeasure;
	}

	/**
	 * @return The tolerance when optimizing for lambda. Used only when lambdaMethod == Opt.
	 */
	double getTol() const {
		return tol;
	}

	/**
	 * Set the tolerance when optimizing for lambda. Used only when lambdaMethod == Opt.
	 * @param tol The new tolerance when optimizing for lambda. Used only when lambdaMethod == Opt.
	 */
	void setTol(double tol) {
		this->tol = tol;
	}

	/**
	 * @return Ratio of negative links used to estimate lambda (used only if lambdaMethod != User).
	 */
	double getLambdaNegEstRatio() const {
		return lambdaNegEstRatio;
	}

	/**
	 * Set the ratio of negative links used to estimate lambda (used only if lambdaMethod != User).
	 * @param lambdaNegEstRatio The new ratio of negative links used to estimate lambda (used only if lambdaMethod != User).
	 */
	void setLambdaNegEstRatio(double lambdaNegEstRatio) {
		this->lambdaNegEstRatio = lambdaNegEstRatio;
	}

	/**
	 * @return Ratio of positive links used to estimate lambda (used only if lambdaMethod != User).
	 */
	double getLambdaPosEstRatio() const {
		return lambdaPosEstRatio;
	}

	/**
	 * Set the ratio of positive links used to estimate lambda (used only if lambdaMethod != User).
	 * @param lambdaPosEstRatio The new ratio of positive links used to estimate lambda (used only if lambdaMethod != User).
	 */
	void setLambdaPosEstRatio(double lambdaPosEstRatio) {
		this->lambdaPosEstRatio = lambdaPosEstRatio;
	}

	/**
	 * @return The method used to compute edge lengths.
	 */
	EdgeLengthMethod getEdgeLengthMethod() const {
		return edgeLengthMethod;
	}

	/**
	 * Set the edge length method.
	 * @param edgeLengthMethod The new edge length method.
	 */
	void setEdgeLengthMethod(EdgeLengthMethod edgeLengthMethod) {
		this->edgeLengthMethod = edgeLengthMethod;
	}

	/**
	 * Destructor.
	 */
	virtual ~UMPSPredictor() = default;

}
;
}
/* namespace LinkPred */

#endif /* UMPSPREDICTOR_HPP_ */
