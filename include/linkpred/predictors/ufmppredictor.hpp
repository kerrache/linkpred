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
 * @brief Contains the implementation of a fast mapping link predictor.
 */

#ifndef UFMPPREDICTOR_HPP_
#define UFMPPREDICTOR_HPP_

#include <linkpred/predictors/ulpredictor.hpp>
#include "linkpred/core/netdistcalculator.hpp"
#include "linkpred/perf/perfmeasure.hpp"
#include "linkpred/utils/log.hpp"
#include <memory>
#include <cmath>
#include <limits>

namespace LinkPred {
/**
 * @brief Fast mapping link predictor.
 * @tparam NetworkT The network type.
 * @tparam EdgesRandomIteratorT A random iterator type used to iterate on edges.
 * @tparam ScoresRandomIteratorT A random iterator type used to iterate on scores.
 */
template<typename NetworkT = UNetwork<>,
		typename EdgesRandomIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::const_iterator,
		typename ScoresRandomIteratorT = typename std::vector<double>::iterator,
		typename EdgesRandomOutputIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::iterator> class UFMPPredictor: public ULPredictor<
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
	 * @brief An enumeration of the different edge length assignment strategies.
	 */
	enum EdgeLengthMethod {
		AUT, /**< Automatic. */
		CST, /**< All edges are assigned the same value, this simply results in a shortest path link predictor. */
		ADA, /**< Adamic Adar strategy. */
		CNE, /**< Common neighbors. */
		HDI, /**< Hub depromoted index strategy. */
		HPI, /**< Hub promoted index strategy. */
		JID, /**< Jackard index. */
		LHN, /**< Leicht-Holme-Newman index strategy. */
		PAT, /**< Preferential attachment strategy. */
		RAL, /**< Resource allocation strategy. */
		SAI, /**< Salton index strategy. */
		SOI /**< Sorensen index strategy. */
	};

	const std::map<std::string, EdgeLengthMethod> strToEdgeLengthMethod = { {
			"AUT", AUT }, { "CST", CST }, { "ADA", ADA }, { "CNE", CNE }, {
			"HDI", HDI }, { "HPI", HPI }, { "JID", JID }, { "LHN", LHN }, {
			"PAT", PAT }, { "RAL", RAL }, { "SAI", SAI }, { "SOI", SOI } }; /**< Mapping of edge length assignment strategies names to their enumeration values. */

	const std::map<EdgeLengthMethod, std::string> edgeLengthMethodToStr = { {
			AUT, "AUT" }, { CST, "CST" }, { ADA, "ADA" }, { CNE, "CNE" }, { HDI,
			"HDI" }, { HPI, "HPI" }, { JID, "JID" }, { LHN, "LHN" }, { PAT,
			"PAT" }, { RAL, "RAL" }, { SAI, "SAI" }, { SOI, "SOI" } }; /**< Mapping of edge length assignment strategies to their names. */

protected:
	long int seed; /**< The random number generator seed. */
	RandomGen rng; /**< The random number generator. */
	Dijkstra<NetworkT, double, std::size_t> dijkstra; /**< Dijkstra's algorithm. */
	typename NetDistCalculator<NetworkT, double, std::size_t>::EdgeLengthMapSP length; /**< The length map. */
	std::shared_ptr<NetDistCalculator<NetworkT, double, std::size_t>> distCalc; /**< The distance calculator. */
	CacheLevel cacheLevel = NodeCache; /**< Cache level for network distances. */
	bool asp = false; /**< Whether to use approximate shortest paths instead of exact ones. */
	double landmarkRatio = 0; /**< The ratio of nodes used as landmarks if asp is on. */
	LandmarkStrategy landmarkStrategy = Random; /**< Landmark positioning strategy if asp is on. */
	std::set<NodeIdType> landmarks; /**< The landmarks. */
	EdgeLengthMethod edgeLengthMethod = AUT; /**< Edge length assignment strategy. */
	double lengthMethodLinksRatio = 0.5; /**< Ratio of links used to select the edge length method (used only if AUT is chosen).*/
	std::shared_ptr<PerfMeasure<>> perfMeasure; /**< The performance measure used to select edge length method (used only if AUT is chosen).*/

	/**
	 * Compute the length of the edge.
	 * @param edge The edge.
	 * @return The edge length using CST strategy.
	 */
	double CSTEdgeLength(EdgeType const & edge) const {
		return 1;
	}

	/**
	 * Compute the length of the edge.
	 * @param edge The edge.
	 * @return The edge length using ADA strategy.
	 */
	double ADAEdgeLength(EdgeType const & edge) const {
		auto srcNode = NetworkT::start(edge);
		auto endNode = NetworkT::end(edge);
		std::vector<NodeIdType> cn;
		net->getCommonNeighbors(srcNode, endNode, std::back_inserter(cn));
		double sum = 0;
		for (auto it = cn.begin(); it != cn.end(); ++it) {
			sum += 1.0 / std::log(net->getDeg(*it));
		}
		return sum + std::numeric_limits<double>::epsilon();
	}

	/**
	 * Compute the length of the edge.
	 * @param edge The edge.
	 * @return The edge length using CNE strategy.
	 */
	double CNEEdgeLength(EdgeType const & edge) const {
		return net->getNbCommonNeighbors(NetworkT::start(edge),
				NetworkT::end(edge)) + std::numeric_limits<double>::epsilon();
	}

	/**
	 * Compute the length of the edge.
	 * @param edge The edge.
	 * @return The edge length using HDI strategy.
	 */
	double HDIEdgeLength(EdgeType const & edge) const {
		auto srcNode = NetworkT::start(edge);
		auto endNode = NetworkT::end(edge);
		return (double) net->getNbCommonNeighbors(srcNode, endNode)
				/ (std::max(net->getDeg(srcNode), net->getDeg(endNode))
						+ std::numeric_limits<double>::epsilon())
				+ std::numeric_limits<double>::epsilon();
	}

	/**
	 * Compute the length of the edge.
	 * @param edge The edge.
	 * @return The edge length using HPI strategy.
	 */
	double HPIEdgeLength(EdgeType const & edge) const {
		auto srcNode = NetworkT::start(edge);
		auto endNode = NetworkT::end(edge);
		return (double) net->getNbCommonNeighbors(srcNode, endNode)
				/ (std::min(net->getDeg(srcNode), net->getDeg(endNode))
						+ std::numeric_limits<double>::epsilon())
				+ std::numeric_limits<double>::epsilon();
	}

	/**
	 * Compute the length of the edge.
	 * @param edge The edge.
	 * @return The edge length using JID strategy.
	 */
	double JIDEdgeLength(EdgeType const & edge) const {
		auto srcNode = NetworkT::start(edge);
		auto endNode = NetworkT::end(edge);
		double degI = net->getDeg(srcNode);
		double degJ = net->getDeg(endNode);
		double nbCN = net->getNbCommonNeighbors(srcNode, endNode);
		return nbCN
				/ (degI + degJ - nbCN + std::numeric_limits<double>::epsilon())
				+ std::numeric_limits<double>::epsilon();

	}

	/**
	 * Compute the length of the edge.
	 * @param edge The edge.
	 * @return The edge length using LHN strategy.
	 */
	double LHNEdgeLength(EdgeType const & edge) const {
		auto srcNode = NetworkT::start(edge);
		auto endNode = NetworkT::end(edge);
		return (double) net->getNbCommonNeighbors(srcNode, endNode)
				/ (net->getDeg(srcNode) * net->getDeg(endNode)
						+ std::numeric_limits<double>::epsilon())
				+ std::numeric_limits<double>::epsilon();
	}

	/**
	 * Compute the length of the edge.
	 * @param edge The edge.
	 * @return The edge length using PAT strategy.
	 */
	double PATEdgeLength(EdgeType const & edge) const {
		double ki = net->getDeg(NetworkT::start(edge));
		double kj = net->getDeg(NetworkT::end(edge));
		return ki * kj + std::numeric_limits<double>::epsilon();
	}

	/**
	 * Compute the length of the edge.
	 * @param edge The edge.
	 * @return The edge length using RAL strategy.
	 */
	double RALEdgeLength(EdgeType const & edge) const {
		auto srcNode = NetworkT::start(edge);
		auto endNode = NetworkT::end(edge);
		std::vector<NodeIdType> cn;
		net->getCommonNeighbors(srcNode, endNode, std::back_inserter(cn));
		double sum = 0;
		for (auto it = cn.begin(); it != cn.end(); ++it) {
			sum += 1.0 / net->getDeg(*it);
		}
		return sum + std::numeric_limits<double>::epsilon();
	}

	/**
	 * Compute the length of the edge.
	 * @param edge The edge.
	 * @return The edge length using SAI strategy.
	 */
	double SAIEdgeLength(EdgeType const & edge) const {
		auto srcNode = NetworkT::start(edge);
		auto endNode = NetworkT::end(edge);
		return (double) net->getNbCommonNeighbors(srcNode, endNode)
				/ (std::sqrt(
						(double) net->getDeg(srcNode) * net->getDeg(endNode))
						+ std::numeric_limits<double>::epsilon())
				+ std::numeric_limits<double>::epsilon();
	}

	/**
	 * Compute the length of the edge.
	 * @param edge The edge.
	 * @return The edge length using SOI strategy.
	 */
	double SOIEdgeLength(EdgeType const & edge) const {
		auto srcNode = NetworkT::start(edge);
		auto endNode = NetworkT::end(edge);
		return 2.0 * net->getNbCommonNeighbors(srcNode, endNode)
				/ (net->getDeg(srcNode) + net->getDeg(endNode)
						+ std::numeric_limits<double>::epsilon())
				+ std::numeric_limits<double>::epsilon();
	}

	/**
	 * Selects automatically the edge length method.
	 */
	void selectEdgeLengthMethod();

	/**
	 * Compute the performance.
	 * @param posBegin Iterator to the first positive link score.
	 * @param posEnd Iterator to one-past-the-last positive link score.
	 * @param negBegin Iterator to the first negative link score.
	 * @param negEnd Iterator to one-past-the-last positive link score.
	 * @return The value of the performance measure.
	 */
	double getPerf(typename std::vector<double>::iterator posBegin,
			typename std::vector<double>::iterator posEnd,
			typename std::vector<double>::iterator negBegin,
			typename std::vector<double>::iterator negEnd) const;

	/**
	 * Create the length map from network topology.
	 */
	void createLengthMap();

	/**
	 * Position the landmarks.
	 */
	void setLandmarks();

public:
	/**
	 * @param net The network.
	 * @param seed The random number generator's seed.
	 */
	UFMPPredictor(std::shared_ptr<NetworkT const> net, long int seed) :
			ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
					EdgesRandomOutputIteratorT>(net), seed(seed), rng(seed), dijkstra(
					net) {
		name = "FMP";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	UFMPPredictor(UFMPPredictor const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	UFMPPredictor & operator =(UFMPPredictor const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	UFMPPredictor(UFMPPredictor && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	UFMPPredictor & operator =(UFMPPredictor && that) = default;

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

		case AUT: {
			throw std::logic_error(
					"Cannot predict with undetermined edge length method");
		}

		case CST: {
			if (net->isEdge(srcNode, endNode)) {
				auto res = distCalc->getIndDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / CSTEdgeLength(e));
			} else {
				auto res = distCalc->getDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / CSTEdgeLength(e));
			}
			break;
		}

		case ADA: {
			if (net->isEdge(srcNode, endNode)) {
				auto res = distCalc->getIndDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / ADAEdgeLength(e));
			} else {
				auto res = distCalc->getDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / ADAEdgeLength(e));
			}
			break;
		}

		case CNE: {
			if (net->isEdge(srcNode, endNode)) {
				auto res = distCalc->getIndDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / CNEEdgeLength(e));
			} else {
				auto res = distCalc->getDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / CNEEdgeLength(e));
			}
			break;
		}

		case HDI: {
			if (net->isEdge(srcNode, endNode)) {
				auto res = distCalc->getIndDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / HDIEdgeLength(e));
			} else {
				auto res = distCalc->getDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / HDIEdgeLength(e));
			}
			break;
		}

		case HPI: {
			if (net->isEdge(srcNode, endNode)) {
				auto res = distCalc->getIndDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / HPIEdgeLength(e));
			} else {
				auto res = distCalc->getDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / HPIEdgeLength(e));
			}
			break;
		}

		case JID: {
			if (net->isEdge(srcNode, endNode)) {
				auto res = distCalc->getIndDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / JIDEdgeLength(e));
			} else {
				auto res = distCalc->getDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / JIDEdgeLength(e));
			}
			break;
		}

		case LHN: {

			if (net->isEdge(srcNode, endNode)) {
				auto res = distCalc->getIndDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / LHNEdgeLength(e));
			} else {
				auto res = distCalc->getDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / LHNEdgeLength(e));
			}
			break;
		}

		case PAT: {

			if (net->isEdge(srcNode, endNode)) {
				auto res = distCalc->getIndDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / PATEdgeLength(e));
			} else {
				auto res = distCalc->getDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / PATEdgeLength(e));
			}
			break;
		}

		case RAL: {

			if (net->isEdge(srcNode, endNode)) {
				auto res = distCalc->getIndDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / RALEdgeLength(e));
			} else {
				auto res = distCalc->getDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / RALEdgeLength(e));
			}
			break;

		}

		case SAI: {

			if (net->isEdge(srcNode, endNode)) {
				auto res = distCalc->getIndDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / SAIEdgeLength(e));
			} else {
				auto res = distCalc->getDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / SAIEdgeLength(e));
			}
			break;
		}

		case SOI: {

			if (net->isEdge(srcNode, endNode)) {
				auto res = distCalc->getIndDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / SOIEdgeLength(e));
			} else {
				auto res = distCalc->getDist(srcNode, endNode);
				sc = 1.0 / (1 + res.first / SOIEdgeLength(e));
			}
			break;
		}

		default:
			throw std::runtime_error("Unknown edge length assignment strategy");
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
	typename NetDistCalculator<NetworkT, double, std::size_t>::EdgeLengthMapSP getLength() const {
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
	 * @return The edge length method.
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
	 * @return Ratio of links used to select the edge length method.
	 */
	double getLengthMethodLinksRatio() const {
		return lengthMethodLinksRatio;
	}

	/**
	 * @return Ratio of links used to select the edge length method.
	 */
	void setLengthMethodLinksRatio(double lengthMethodLinksRatio) {
		if (lengthMethodLinksRatio <= 0 || lengthMethodLinksRatio > 1) {
			throw std::invalid_argument(
					"The ratio of links used to select the edge length method must be between 0 exclusive, and 1 inclusive");
		}
		this->lengthMethodLinksRatio = lengthMethodLinksRatio;
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
	 * Destructor.
	 */
	virtual ~UFMPPredictor() = default;

};

}
/* namespace LinkPred */

#endif /* UFMPPREDICTOR_HPP_ */
