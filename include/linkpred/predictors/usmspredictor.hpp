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
 * @brief Contains the implementation of a single MDS run link predictor.
 */

#ifndef USMSPREDICTOR_HPP_
#define USMSPREDICTOR_HPP_

#include <linkpred/predictors/ulpredictor.hpp>
#include <cmath>

namespace LinkPred {
/**
 * @brief Link predictor that uses a single MDS run (with constant degrees).
 * @tparam NetworkT The network type.
 * @tparam EdgesRandomIteratorT A random iterator type used to iterate on edges.
 * @tparam ScoresRandomIteratorT A random iterator type used to iterate on scores.
 */
template<typename NetworkT = UNetwork<>,
		typename EdgesRandomIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::const_iterator,
		typename ScoresRandomIteratorT = typename std::vector<double>::iterator,
		typename EdgesRandomOutputIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::iterator> class USMSPredictor: public ULPredictor<
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
	 * Enumeration of the possible link statuses
	 */
	enum LinkStatus {
		Pos, /**< Positive link. */
		Neg, /**< Negative link. */
		Unk, /**< Unknown. */
	};

protected:
	long int seed; /**< Random generator seed. */
	RandomGen rng; /**< Random number generator. */
	double hProb = 0.95; /**< High probability. This is assigned to connected couples. */
	double lProb = 0.05; /**< Low probability. This is assigned to disconnected couples. */
	double alpha = 2.5; /**< The parameter alpha. */
	std::size_t dim = 2; /**< The dimension. */
	std::size_t nbNodes = 0; /**< Number of nodes. */
	std::size_t nbCouples = 0; /**< Number of couples. */
	std::vector<double> coords; /**< The coordinates of the nodes. */
	std::vector<double> deg; /**< The (expected) degrees of the nodes. */
	std::vector<double> degProd; /**< Product of the degrees. */
	std::vector<double> sqDist; /**< Squared distances between nodes. */
	std::vector<LinkStatus> linkStatus; /**< Link statuses. */
	std::vector<double> weight; /**< Link weight in MDS. */
	double wPos = 1; /**< Weight of positive links. */
	double wNeg = 1; /**< Weight of negative links. */
	std::size_t nbMdsRuns = 5; /**< Number of times MDS is run. */

	/**
	 * Compute the distances.
	 */
	void compDist();

	/**
	 * Solve the MDS problem.
	 */
	void mds();

public:

	/**
	 * @param net The network.
	 * @param seed The random number generator's seed.
	 */
	USMSPredictor(std::shared_ptr<NetworkT const> net, long int seed) :
			ULPredictor<NetworkT, EdgesRandomIteratorT>(net), seed(seed), rng(
					seed) {
		name = "SMS";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	USMSPredictor(USMSPredictor const & that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	USMSPredictor & operator =(USMSPredictor const & that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	USMSPredictor(USMSPredictor && that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	USMSPredictor & operator =(USMSPredictor && that) = default;

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
		return std::pow(1 + std::sqrt(sqDist[k]) / degProd[k], -alpha);
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
	 * @return alpha.
	 */
	double getAlpha() const {
		return alpha;
	}

	/**
	 * @param alpha The new value of alpha.
	 */
	void setAlpha(double alpha) {
		this->alpha = alpha;
	}

	/**
	 * @return hProb.
	 */
	double getHProb() const {
		return hProb;
	}

	/**
	 * @param hProb The new value of hProb.
	 */
	void setHProb(double hProb) {
		this->hProb = hProb;
	}

	/**
	 * @return lProb.
	 */
	double getLProb() const {
		return lProb;
	}

	/**
	 * @param lProb The new value of lProb.
	 */
	void setLProb(double lProb) {
		this->lProb = lProb;
	}

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
	 * Destructor.
	 */
	virtual ~USMSPredictor() = default;

};

}
/* namespace LinkPred */

#endif /* USMSPREDICTOR_HPP_ */
