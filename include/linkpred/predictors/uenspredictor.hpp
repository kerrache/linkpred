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
 * @brief Contains the implementation of an ensemble link predictor.
 */

#ifndef UENSPREDICTOR_HPP_
#define UENSPREDICTOR_HPP_

#include <linkpred/predictors/ulpredictor.hpp>
#include "linkpred/perf/perfmeasure.hpp"
#include "linkpred/utils/log.hpp"
#include <memory>
#include <cmath>
#include <limits>

namespace LinkPred {
// TODO check literature on ensemble learning and transform this class into a generic ensemble predictor
/**
 * @brief Ensemble link predictor.
 * @tparam NetworkT The network type.
 * @tparam EdgesRandomIteratorT A random iterator type used to iterate on edges.
 * @tparam ScoresRandomIteratorT A random iterator type used to iterate on scores.
 */
template<typename NetworkT = UNetwork<>,
		typename EdgesRandomIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::const_iterator,
		typename ScoresRandomIteratorT = typename std::vector<double>::iterator,
		typename EdgesRandomOutputIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::iterator> class UENSPredictor: public ULPredictor<
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
	 * @brief An enumeration of the different edge length assignment strategies.
	 */
	enum EdgeScoreMethod {
		AUT, /**< Automatic. */
		CST, /**< All edges are assigned the same value. */
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

	const std::map<std::string, EdgeScoreMethod> strToEdgeScoreMethod = { {
			"AUT", AUT }, { "CST", CST }, { "ADA", ADA }, { "CNE", CNE }, {
			"HDI", HDI }, { "HPI", HPI }, { "JID", JID }, { "LHN", LHN }, {
			"PAT", PAT }, { "RAL", RAL }, { "SAI", SAI }, { "SOI", SOI } }; /**< Mapping of edge length assignment strategies names to their enumeration values. */

	const std::map<EdgeScoreMethod, std::string> edgeScoreMethodToStr = { { AUT,
			"AUT" }, { CST, "CST" }, { ADA, "ADA" }, { CNE, "CNE" }, { HDI,
			"HDI" }, { HPI, "HPI" }, { JID, "JID" }, { LHN, "LHN" }, { PAT,
			"PAT" }, { RAL, "RAL" }, { SAI, "SAI" }, { SOI, "SOI" } }; /**< Mapping of edge length assignment strategies to their names. */

protected:
	long int seed = 0; /**< The random number generator seed. */
	RandomGen rng; /**< The random number generator. */
	EdgeScoreMethod scoreMethod = AUT; /**< Edge length assignment strategy. */
	double scoreMethodLinksRatio = 0.5; /**< Ratio of links used to select the edge length method (used only if AUT is chosen).*/
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
	void selectEdgeScoreMethod();

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

public:
	/**
	 * @param net The network.
	 * @param seed The random number generator's seed.
	 */
	UENSPredictor(std::shared_ptr<NetworkT const> net, long int seed) :
			ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
					EdgesRandomOutputIteratorT>(net), seed(seed), rng(seed) {
		name = "ENS";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	UENSPredictor(UENSPredictor const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	UENSPredictor & operator =(UENSPredictor const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	UENSPredictor(UENSPredictor && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	UENSPredictor & operator =(UENSPredictor && that) = default;

	/**
	 * Initialize the predictor.
	 */
	virtual void init();

	/**
	 * Learn.
	 */
	virtual void learn();

	/**
	 * Predict the links.
	 * @param begin Beginning of the links to be predicted.
	 * @param end end of the links to be predicted.
	 * @param scores Beginning of scores.
	 */
	virtual void predict(EdgesRandomIteratorT begin, EdgesRandomIteratorT end,
			ScoresRandomIteratorT scores);

	/**
	 * @return The edge score method.
	 */
	EdgeScoreMethod getScoreMethod() const {
		return scoreMethod;
	}

	/**
	 * Set the edge score method.
	 * @param edgeScoreMethod The new edge score method.
	 */
	void setScoreMethod(EdgeScoreMethod edgeScoreMethod) {
		this->scoreMethod = edgeScoreMethod;
	}

	/**
	 * @return The ratio of links used to determine the best edge score method.
	 */
	double getScoreMethodLinksRatio() const {
		return scoreMethodLinksRatio;
	}

	/**
	 * Set the ratio of links used to determine the best edge score method.
	 * @param scoreMethodLinksRatio The new ratio of links used to determine the best edge score method.
	 */
	void setScoreMethodLinksRatio(double scoreMethodLinksRatio) {
		if (scoreMethodLinksRatio <= 0 || scoreMethodLinksRatio > 1) {
			throw std::invalid_argument(
					"The ratio of links used to select the edge length method must be between 0 exclusive, and 1 inclusive");
		}
		this->scoreMethodLinksRatio = scoreMethodLinksRatio;
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
	 * Destructor.
	 */
	virtual ~UENSPredictor() = default;

};

}
/* namespace LinkPred */

#endif /* UENSPREDICTOR_HPP_ */
