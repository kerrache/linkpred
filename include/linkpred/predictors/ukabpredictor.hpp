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
 * @brief Contains the implementation of a scalable popularity-similarity link predictor.
 */

#ifndef UKABPREDICTOR_HPP_
#define UKABPREDICTOR_HPP_

#include <linkpred/predictors/ukabpredictor/kablambdacg.hpp>
#include <linkpred/predictors/ulpredictor.hpp>
#include "linkpred/core/netdistcalculator.hpp"
#include "linkpred/perf/perfmeasure.hpp"
#include "linkpred/utils/log.hpp"
#include <memory>
#include <cmath>
#include <limits>
#include <iostream>

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
				typename NetworkT::EdgeType>::iterator> class UKABPredictor: public ULPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT> {

public:
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::net;
#ifdef WITH_OPENMP
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
	EdgesRandomOutputIteratorT>::parallel; /**< Whether the predictor runs in parallel. */
#endif
#ifdef WITH_MPI
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
	EdgesRandomOutputIteratorT>::comm; /**< The MPI communicator. */
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
	EdgesRandomOutputIteratorT>::distributed; /**< Enable/disable distributed parallelism. */
#endif
	using ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
			EdgesRandomOutputIteratorT>::name;
	using NodeIdType = typename ULPredictor<NetworkT, EdgesRandomIteratorT,ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::NodeIdType;
	using EdgeType = typename ULPredictor<NetworkT, EdgesRandomIteratorT,ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::EdgeType;

	/**
	 * Lambda method.
	 */
	enum LambdaMethT {
		User, Est, Opt, Scan
	};

	/**
	 * Edge length type.
	 */
	enum EdgeLengthT {
		One, Sum, Prod, NSC, SNSC, RA1, RA2

	};

	/**
	 * Score type.
	 */
	enum ScoreT {
		Dist, SumDist, ProdDist, SNSCSc, SNSCDist, RA1Dist, RA2Dist
	};

protected:

	long int seed = 0; /**< The random number generator seed. */
	RandomGen rng; /**< The random number generator. */
	Dijkstra<NetworkT, double, std::size_t> dijkstra; /**< Dijkstra's algorithm. */
	typename ESPLDistCalculator<NetworkT, double, std::size_t>::EdgeLengthMapSP length; /**< The length map. */
	std::shared_ptr<ESPLDistCalculator<NetworkT, double, std::size_t>> distCalc; /**< The distance calculator. */
	CacheLevel cacheLevel = NetworkCache; /**< The cache level. */
	std::size_t maxDeg = 0; /**< The maximum degree in the network. */
	std::size_t lim = 2; /**< Limit on the number of hops. */
	bool useHops = true; /**< Use hops when computing similarity. If false, distance is used instead. */
	double lambda = 0.5;
	std::vector<double> lambdas =
			{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
	LambdaMethT lambdaMeth = Est;
	std::size_t sampleSize = 10000; /**< Sample size for finding lambda. */
	EdgeLengthT edgeLengthT = SNSC;
	ScoreT scoreT = SNSCDist;
	double zeta = -1;
	inline double phi(double k) const {
		return std::log(k + 1);
	}

	inline double sum(double ki, double kj) const {
//		return (phi(ki) + phi(kj)) / (2 * phi(maxDeg));
		return (phi(ki) + phi(kj)) / (2 * phi(maxDeg));
	}

	inline double prod(double ki, double kj) const {
		return std::sqrt((phi(ki) / phi(maxDeg)) * (phi(kj) / phi(maxDeg)));
	}

	inline double nsc(NodeIdType srcNode, NodeIdType endNode, double ki,
			double kj) const {
		std::vector<NodeIdType> cn;
		net->getCommonNeighbors(srcNode, endNode, std::back_inserter(cn));
		double ns = 1;
		for (auto it = cn.begin(); it != cn.end(); ++it) {
			double kk = net->getDeg(*it);
			ns *= phi(kk) / phi(maxDeg);
		}
		return 1 - ns;
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
	 * Compute the length of the edge.
	 * @param edge The edge.
	 * @return The edge length using PAT method.
	 */
	inline double edgeLength(EdgeType const & edge) const {
		auto i = NetworkT::start(edge);
		auto j = NetworkT::end(edge);

		double ki = net->getDeg(i);
		double kj = net->getDeg(j);
		switch (edgeLengthT) {
		case One:
			return 1.0;
		case Sum:
			return sum(ki, kj);
		case Prod:
			return prod(ki, kj);
		case NSC:
			return nsc(i, j, ki, kj);
		case SNSC: {
			double ns = nsc(i, j, ki, kj);
//			return sum(ki, kj);
//			return 1 / (1+sum(ki, kj) + ns); //l=(2.07,6), s=(1.88,5)
//			return sum(ki, kj) / (1 - zeta + ns * ns); //l=(2.07,6), s=(1.88,5)
//			return sum(ki, kj) / (1 - zeta + (zeta) * ns); //l=(2.000, 6), s=(1.73, 7)
//			std::cout << net->getLabel(i) << "\t" << net->getLabel(j) << "\t"
//					<< "sum: " << sum(ki, kj) << " nsc: " << ns << std::endl;
			return sum(ki, kj) / (zeta + (1 - zeta) * ns); //l=(1.9706, 7), s=(1.9048, 7)
//			return sum(ki, kj) / (lambda + (1 - lambda) * ns); //l=(2.0476,7), s=(1.7647,7)
//			return sum(ki, kj) / (1 + ns); //l=(2.0476,7), s=(1.7647,7)
//			return sum(ki, kj) / (1 + ns); //l=(2.0476,7), s=(1.7647,7)
//			return std::log2(1 + sum(ki, kj)) / (0.5 + std::log2(1 + ns)); //l=(2.0952, 6), s=(2.1176, 4)
//			return std::log2(1 + sum(ki, kj)) / (1 + std::log2(1 + ns)); //l=(2.2143, 7), s=(1.8529, 6)
//			return std::log2(1 + sum(ki, kj)) / (lambda + std::log2(1 + ns)); //lambda=formula, l=(1.9762, 8), s=(2.1765, 3)
//			return std::log2(1 + sum(ki, kj)) / (zeta + std::log2(1 + ns)); //zeta=formula, l=(2.0000, 9), s=(2.0588, 4)
//			return std::pow(sum(ki, kj), zeta + ns); //zeta =cc, l=(2.1905, 8), s=(2.0588, 5)
//			return std::pow(sum(ki, kj), ns); //l=(2.0476, 7), s=(2.0294, 5)
//			return std::log2(1 + sum(ki, kj)) / std::log2(1 + ns);
		}
//			return sum(ki, kj) / (1 + nsc(srcNode, endNode, ki, kj));
//			return sum(ki, kj)
//					/ (lambda + (1 - lambda) * nsc(srcNode, endNode, ki, kj));
//			return sum(ki, kj)
//					/ (1 + (1 - lambda) * nsc(srcNode, endNode, ki, kj));

//			return std::log2(1+sum(ki, kj)) / std::log2(1 + nsc(srcNode, endNode, ki, kj));
//			return sum(ki, kj)
//					/ (lambda / (1 - lambda) + nsc(srcNode, endNode, ki, kj));
//			return std::pow(sum(ki, kj), lambda)
//					* std::pow(nsc(srcNode, endNode, ki, kj), 1 - lambda);
//			return sum(ki, kj) / (1.0 + nsc(srcNode, endNode, ki, kj));
//			return 1.0 + (1 - lambda) * sum(ki, kj)
//					- lambda * nsc(srcNode, endNode, ki, kj);
//			return (1 - lambda) * sum(ki, kj)
//					/ (1 + lambda * nsc(srcNode, endNode, ki, kj));
//			return lambda / (1 - lambda) * sum(ki, kj) / nsc(srcNode, endNode, ki, kj);
//			return ((1 - lambda) * sum(ki, kj) + lambda * nsc(srcNode, endNode, ki, kj));

		case RA1:
			return RA1EdgeLength(edge);

		case RA2:
			return RA2EdgeLength(edge);

		default:
			throw std::runtime_error("Unknown edge length type");
		}
	}

	/**
	 * Create the length map from network topology.
	 */
	void createLengthMap();

	/**
	 * @return A sparse node map containing the distance to all nodes having a nonzero similarity with srcNode.
	 * @param srcNode Source node.
	 */
	typename NetworkT::template NodeSMapSP<std::pair<double, size_t>> getFiniteDistMap(
			NodeIdType srcNode) {
		auto distMap = distCalc->getFinDistMapNoNeighb(srcNode);
		for (auto it = distMap->begin(); it != distMap->end(); ++it) {
			it->second.first = 1.0 + it->second.first;
		}
		return distMap;
	}

	/**
	 * Computes the similarity between two nodes.
	 * @param srcNode Source node.
	 * @param endNode End node.
	 */
	double getDist(NodeIdType srcNode, NodeIdType endNode) {
		double dist = 1 + distCalc->getDist(srcNode, endNode).first;
		return dist;
	}

	/**
	 * Computes the similarity between two nodes.
	 * @param srcNode Source node.
	 * @param endNode End node.
	 */
	double getDist(EdgeType const & e) {
		return getDist(net->start(e), net->end(e));
	}

	/**
	 * Compute the score of a single edge.
	 * @param e The edge.
	 * @param dist Distance associated with e.
	 * @return The score of e.
	 */
	double score(EdgeType const & e, double dist);

	/**
	 * Find lambda.
	 */
	void findLambda();

	/**
	 * Scan for the best value of lambda.
	 */
	void lambdaScan();

	/**
	 * Find lambda using optimization.
	 */
	void lambdaOpt();

public:
	/**
	 * @param net The network.
	 * @param seed The random number generator's seed.
	 */
	UKABPredictor(std::shared_ptr<NetworkT const> net, long int seed) :
			ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
					EdgesRandomOutputIteratorT>(net), seed(seed), rng(seed), dijkstra(
					net) {
		name = "KAB";
		maxDeg = net->getMaxDeg();
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	UKABPredictor(UKABPredictor const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	UKABPredictor & operator =(UKABPredictor const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	UKABPredictor(UKABPredictor && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	UKABPredictor & operator =(UKABPredictor && that) = default;

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
	virtual double score(EdgeType const & e);

	/**
	 * Predict the links.
	 * @param begin Beginning of the links to be predicted.
	 * @param end end of the links to be predicted.
	 * @param scores Beginning of scores.
	 */
	virtual void predict(EdgesRandomIteratorT begin, EdgesRandomIteratorT end,
			ScoresRandomIteratorT scores);

	/**
	 * Finds the k negative edges with the top score. Ties are broken randomly.
	 * @param k The number of edges to find.
	 * @param eit An output iterator where the edges are written.
	 * @param sit An output iterator where the scores are written. The scores are written in the same order as the edges.
	 * @return The number of negative edges inserted. It is the minimum between k and the number of negative edges in the network.
	 */
	virtual std::size_t top(std::size_t k, EdgesRandomOutputIteratorT eit,
			ScoresRandomIteratorT sit);

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

	std::size_t getLim() const {
		return lim;
	}

	void setLim(std::size_t lim) {
		this->lim = lim;
	}

	/**
	 * @return True if the number of hops is used to compute similarity.
	 */
	bool isUseHops() const {
		return useHops;
	}

	/**
	 * @param useHops Set whether the number of hops is used to compute similarity.
	 */
	void setUseHops(bool useHops) {
		this->useHops = useHops;
	}

	double getLambda() const {
		return lambda;
	}

	void setLambda(double lambda = 0.25) {
		this->lambda = lambda;
	}

	EdgeLengthT getEdgeLengthT() const {
		return edgeLengthT;
	}

	void setEdgeLengthT(EdgeLengthT edgeLengthT) {
		this->edgeLengthT = edgeLengthT;
	}

	ScoreT getScoreT() const {
		return scoreT;
	}

	void setScoreT(ScoreT scoreT) {
		this->scoreT = scoreT;
	}

	const std::vector<double>& getLambdas() const {
		return lambdas;
	}

	void setLambdas(const std::vector<double>& lambdas) {
		this->lambdas = lambdas;
	}

	LambdaMethT getLambdaMeth() const {
		return lambdaMeth;
	}

	void setLambdaMeth(LambdaMethT lambdaMeth) {
		this->lambdaMeth = lambdaMeth;
	}

	double getZeta() const {
		return zeta;
	}

	void setZeta(double zeta = 1) {
		this->zeta = zeta;
	}

	/**
	 * Destructor.
	 */
	virtual ~UKABPredictor() = default;

};
}
/* namespace LinkPred */

#endif /* UKABPREDICTOR_HPP_ */
