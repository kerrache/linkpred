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
 * @ingroup Predictors
 * @brief Contains the implementation of the hypermap link predictor.
 */

#ifndef UHYPPREDICTOR_HPP_
#define UHYPPREDICTOR_HPP_

#include <linkpred/predictors/undirected/ulpredictor.hpp>
#include "linkpred/utils/randomgen.hpp"
#include <vector>
#include <memory>

namespace LinkPred {
constexpr double MathPI = 3.141592653589793238462643383279502884L; /**< Pi. */

/**
 * @brief Hypermap predictor.
 * @details This is a modified and wrapped version of the code provided by the authors.
 * @tparam Network The network type.
 * @tparam EdgeRndIt A random iterator type used to iterate on edges.
 * @tparam ScoreRndIt A random iterator type used to iterate on scores.
 */
template<typename Network = UNetwork<>,
		typename EdgeRndIt = typename std::vector<typename Network::Edge>::const_iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator,
		typename EdgeRndOutIt = typename std::vector<typename Network::Edge>::iterator> class UHYPPredictor: public ULPredictor<
		Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt> {

	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::net; /**< The network. */
#ifdef LINKPRED_WITH_OPENMP
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::parallel; /**< Whether the predictor runs in parallel. */
#endif
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::name; /**< The name of the predictor. */
	using NodeID = typename ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::NodeID; /**< The node IDs type. */
	using Edge = typename ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::Edge; /**< The edges type. */

	/**
	 * @brief Graph node.
	 */
	class Node {

	private:
		NodeID id; /**< Node ID. */
		int flag; /**< Flag. */
		double theta = 0; /**< The angle coordinate. */
		double radius = 0; /**< The distance from the center. */
		double radius_init = 0; /**< Initial radius. */
		double R_i = 0; /**< R_i. */
		std::vector<std::shared_ptr<Node>> adjNodeList; /**< Adjacent nodes. */

	public:
		/**
		 * Constructor.
		 * @param id The node ID.
		 */
		Node(NodeID const & id) {
			this->id = id;
			flag = 0;
		}

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		Node(Node const & that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		Node & operator =(Node const & that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		Node(Node && that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		Node & operator =(Node && that) = default;

		/**
		 * Set The angle coordinate.
		 * @param angle The new angle coordinate.
		 */
		void setAngle(double angle) {
			theta = angle;
		}

		/**
		 * @return The angle coordinate.
		 */
		double getAngle() {
			return theta;
		}

		/**
		 * Set the radius coordinate.
		 * @param r The new radius coordinate.
		 */
		void setRadius(double r) {
			radius = r;
		}

		/**
		 * @return The angle coordinate.
		 */
		double getRadius() {
			return radius;
		}

		/**
		 * Set the initial radius.
		 * @param r The new initial radius.
		 */
		void setInitRadius(double r) {
			radius_init = r;
		}

		/**
		 * @return The initial radius.
		 */
		double getInitRadius() {
			return radius_init;
		}

		/**
		 * Set R.
		 * @param R The new value of R.
		 */
		void set_R(double R) {
			R_i = R;
		}

		/**
		 * @return R.
		 */
		double get_R() {
			return R_i;
		}

		/**
		 * Set the node flag.
		 * @param f The new value of flag.
		 */
		void setFlag(int f) {
			flag = f;
		}

		/**
		 * @return The node flag.
		 */
		int getFlag() {
			return flag;
		}

		/**
		 * @return The node ID.
		 */
		NodeID const & getId() {
			return id;
		}

		/**
		 * Add an adjacent node.
		 * @param adj A pointer to the node.
		 */
		void addAdjNode(std::shared_ptr<Node> adj) {
			adjNodeList.push_back(adj);
		}

		/**
		 * @return The list of adjacent nodes.
		 */
		std::vector<std::shared_ptr<Node>> getAdjNodeList() {
			return adjNodeList;
		}

		/**
		 * Check if a node is adjacent to this.
		 * @param x The node that should be checked.
		 * @return True if x is adjacent to the current node.
		 */
		bool checkLink(std::shared_ptr<Node> x) {
			for (std::size_t i = 0; i < adjNodeList.size(); i++) {
				if (adjNodeList[i] == x) {
					return true;
				}
			}
			return false;
		}

		/**
		 * Destructor.
		 */
		virtual ~Node() = default;
	};

	/**
	 * @brief This class represents a graph.
	 */
	class Graph {

	private:
		std::vector<std::shared_ptr<Node>> nodeList; /**< The list of nodes. */

	public:

		/**
		 * Constructor.
		 */
		Graph() = default;

		/**
		 * Copy constructor.
		 * @param that The object to copy.
		 */
		Graph(Graph const & that) = default;

		/**
		 * Copy assignment operator.
		 * @param that The object to copy.
		 */
		Graph & operator =(Graph const & that) = default;

		/**
		 * Move constructor.
		 * @param that The object to move.
		 */
		Graph(Graph && that) = default;

		/**
		 * Move assignment operator.
		 * @param that The object to move.
		 */
		Graph & operator =(Graph && that) = default;

		/**
		 * Add a node to the graph.
		 * @param nNode The new node.
		 */
		void addNewNode(std::shared_ptr<Node> nNode) {
			nodeList.push_back(nNode);
		}

		/**
		 * Lookup node by ID.
		 * @param id The node ID.
		 * @return The node object having id as its ID.
		 */
		std::shared_ptr<Node> findNode(NodeID id) {
			for (std::size_t i = 0; i < nodeList.size(); i++) {
				if (nodeList[i]->getId() == id)
					return nodeList[i];
			}
			return nullptr;
		}

		/**
		 * @return The node list.
		 */
		std::vector<std::shared_ptr<Node>> getNodeList() {
			return nodeList;
		}

		/**
		 * Destructor.
		 */
		virtual ~Graph() = default;
	};

protected:
	long int seed; /**< The random number generator seed. */
	RandomGen rng; /**< The random number generator. */
	std::map<Edge, double> scores; /**< Links scores. */
	std::shared_ptr<Graph> G; /**< The graph. */
	std::size_t maxDeg = 0; /**< Maximum degree. */
	//Change this according to parameters of network to embed.
	double m = 1.5; /**< The parameter m (see the algorithm description). */
	double L = 1; /**< The parameter L (see the algorithm description). */
	double gamma = 2.1; /**< The power law exponent gamma (see the algorithm description). */
	double zeta = 1; /**< The parameter zeta (see the algorithm description). */
	double T = 0.8; /**< The parameter T (see the algorithm description). */
	std::map<NodeID, std::pair<double, double>> nodesCoord; /**< Nodes coordinates. */

public:

	/**
	 * @param net The network.
	 * @param seed The random number generator's seed.
	 */
	UHYPPredictor(std::shared_ptr<Network const> net, long int seed) :
			ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>(net), seed(
					seed), rng(seed) {
		name = "HYP";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	UHYPPredictor(UHYPPredictor const & that) = delete;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	UHYPPredictor & operator =(UHYPPredictor const & that) = delete;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	UHYPPredictor(UHYPPredictor && that) = delete;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	UHYPPredictor & operator =(UHYPPredictor && that) = delete;

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
	virtual double score(Edge const & e) {
		auto coordsi = nodesCoord.at(Network::start(e));
		auto coordsj = nodesCoord.at(Network::end(e));

		double dtheta = std::abs(coordsi.first - coordsj.first);
		if (dtheta > MathPI) {
			dtheta = 2 * MathPI - dtheta;
		}
		double ri = coordsi.second;
		double rj = coordsj.second;
		double dist = (1 / zeta)
				* std::acosh(
						(std::cosh(zeta * ri) * std::cosh(zeta * rj))
								- (std::sinh(zeta * ri) * std::sinh(zeta * rj)
										* std::cos(dtheta)));
		return 1 / (1 + std::exp((zeta / (2 * T)) * dist));
	}

	/**
	 * Predict the links.
	 * @param begin Beginning of the links to be predicted.
	 * @param end end of the links to be predicted.
	 * @param scores Beginning of scores.
	 */
	virtual void predict(EdgeRndIt begin, EdgeRndIt end, ScoreRndIt scores);

	/**
	 * @return The random number generator seed.
	 */
	long int getSeed() const {
		return seed;
	}

	/**
	 * @return The power law exponent gamma (see the algorithm description).
	 */
	double getGamma() const {
		return gamma;
	}

	/**
	 * Set the power law exponent gamma (see the algorithm description).
	 * @param gamma The new power law exponent gamma (see the algorithm description).
	 */
	void setGamma(double gamma) {
		this->gamma = gamma;
	}

	/**
	 * @return The parameter L (see the algorithm description).
	 */
	double getL() const {
		return L;
	}

	/**
	 * Set the parameter L (see the algorithm description).
	 * @param L The new value of the parameter L (see the algorithm description).
	 */
	void setL(double L) {
		this->L = L;
	}

	/**
	 * @return The parameter m (see the algorithm description).
	 */
	double getM() const {
		return m;
	}

	/**
	 * Set the parameter m (see the algorithm description).
	 * @param m The new value of the parameter m (see the algorithm description).
	 */
	void setM(double m) {
		this->m = m;
	}

	/**
	 * @return The parameter L (see the algorithm description).
	 */
	double getT() const {
		return T;
	}

	/**
	 * Set the parameter L (see the algorithm description).
	 * @param T The new value of the parameter L (see the algorithm description).
	 */
	void setT(double T) {
		this->T = T;
	}

	/**
	 * @return The parameter zeta (see the algorithm description).
	 */
	double getZeta() const {
		return zeta;
	}

	/**
	 * Set the parameter zeta (see the algorithm description).
	 * @param zeta The new value of the parameter zeta (see the algorithm description).
	 */
	void setZeta(double zeta) {
		this->zeta = zeta;
	}

	/**
	 * Destructor.
	 */
	virtual ~UHYPPredictor() = default;
};

} /* namespace LinkPred */

#endif /* UHYPPREDICTOR_HPP_ */
