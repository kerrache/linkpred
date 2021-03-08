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
 * @ingroup GraphAlg
 * @brief Contains the implementation of an algorithm for embedding a network
 * using a a hidden metric space model.
 */

#ifndef HMSM_HPP_
#define HMSM_HPP_

#include "linkpred/graphalg/encoders/encoder.hpp"
#include "linkpred/utils/randomgen.hpp"
#include <cmath>

namespace LinkPred {

/**
 * @brief Contains the implementation of an algorithm for embedding a network
 * using a a hidden metric space model.
 * Reference: R. Alharbi, H. Benhidour, and S. Kerrache. “Link Prediction in
 * Complex Net-works Based on a Hidden Variables Model”. In: 2016 UKSim-AMSS
 * 18th Inter-national Conference on Computer Modelling and Simulation (UKSim).
 * 2016,pages 119–124
 * @tparam Network The network type.
 */
template<typename Network = UNetwork<>> class HMSM: public Encoder<Network> {
public:
	using typename Encoder<Network>::NodeID; /**< Node ID type. */
	using typename Encoder<Network>::CodeMap; /**< Code map type. */
	using typename Encoder<Network>::CodeMapSP; /**< Shared pointer to a node code map. */
	using typename Encoder<Network>::WeightMap; /**< Edge weight map type. */
	using typename Encoder<Network>::WeightMapSP; /**< Shared pointer to an edge weight map. */

protected:
	using Encoder<Network>::net; /**< The network. */
	using Encoder<Network>::name; /**< The name of the encoder. */
	using Encoder<Network>::codeMap; /**< Code map. */
	using Encoder<Network>::dim; /**< Dimension of the embedding. */
#ifdef LINKPRED_WITH_OPENMP
	using Encoder<Network>::parallel; /**< Enable/disable shared-memory parallelism. */
#endif
#ifdef LINKPRED_WITH_MPI
	using Encoder<Network>::distributed; /**< Enable/disable distributed parallelism. */
	using Encoder<Network>::comm; /**< The MPI communicator. */
#endif

	double hProb = 0.95; /**< High probability. This is assigned to connected couples. */
	double lProb = 0.05; /**< Low probability. This is assigned to disconnected couples. */
	double alpha = 2.5; /**< The parameter alpha. */
	int nbMdsRuns = 1; /**< Number of times MDS is run. */
	double tol = 1.0e-5; /**< Tolerance for stopping MDS. */
	double eps = 0.001; /**< Degree product when of the nodes has degree 0. */
	const int DefaultDim = 3; /**< Default embedding dimension. */

private:
	RandomGen rng; /**< Random number generator. */
	int eucDim = 2; /**< Dimension of the Euclidean space = dim -1. */

public:
	/**
	 * Constructor.
	 * @param net The network.
	 * @param seed Random number generator seed.
	 */
	HMSM(std::shared_ptr<Network const> net, long int seed) :
			Encoder<Network>(net), rng(seed) {
		dim = DefaultDim;
		name = "HMS";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	HMSM(HMSM const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	HMSM& operator =(HMSM const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	HMSM(HMSM &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	HMSM& operator =(HMSM &&that) = default;

	/**
	 * Initialize encoder.
	 */
	virtual void init();

	/**
	 * Encode the network.
	 */
	virtual void encode();

	/**
	 * Set edge weight map.
	 * @param weightMap Edge weight map.
	 */
	virtual void setWeightMap(const WeightMapSP &weightMap) {
		throw std::runtime_error(
				"HMSM encoder does not support weighted networks.");
	}

	/**
	 * @return The default embedding dimension.
	 */
	const int getDefaultDim() const {
		return DefaultDim;
	}

	/**
	 * @return The parameter alpha of the HMSM model.
	 */
	double getAlpha() const {
		return alpha;
	}

	/**
	 * Set the parameter alpha of the HMSM model.
	 * @param alpha The parameter alpha of the HMSM model.
	 */
	void setAlpha(double alpha) {
		this->alpha = alpha;
	}

	/**
	 * @return Degree product when of the nodes has degree 0.
	 */
	double getEps() const {
		return eps;
	}

	/**
	 * Set degree product when of the nodes has degree 0.
	 * @param eps Degree product when of the nodes has degree 0.
	 */
	void setEps(double eps) {
		this->eps = eps;
	}

	/**
	 * Probability assigned to connected couples.
	 */
	double getHProb() const {
		return hProb;
	}

	/**
	 * Set the probability assigned to connected couples.
	 * @param hProb Probability assigned to connected couples.
	 */
	void setHProb(double hProb) {
		this->hProb = hProb;
	}

	/**
	 * Probability assigned to disconnected couples.
	 */
	double getLProb() const {
		return lProb;
	}

	/**
	 * Set the probability assigned to disconnected couples.
	 * @param lProb Probability assigned to disconnected couples.
	 */
	void setProb(double lProb) {
		this->lProb = lProb;
	}

	/**
	 * @return The number of times MDS is run (with default starting point).
	 */
	int getNbMdsRuns() const {
		return nbMdsRuns;
	}

	/**
	 * Set the number of times MDS is run (with default starting point).
	 * @param nbMdsRuns Number of times MDS is run (with default starting point).
	 */
	void setNbMdsRuns(int nbMdsRuns) {
		this->nbMdsRuns = nbMdsRuns;
	}

	/**
	 * @return Tolerance for stopping the MDS optimization problem.
	 */
	double getTol() const {
		return tol;
	}

	/**
	 * Set the tolerance for stopping the MDS optimization problem.
	 * @param tol Tolerance for stopping the MDS optimization problem.
	 */
	void setTol(double tol) {
		this->tol = tol;
	}

	/**
	 * Destructor.
	 */
	virtual ~HMSM() = default;

};

} /* namespace LinkPred */

#endif /* HMSM_HPP_ */
