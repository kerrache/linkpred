/*
 * This file is part of LinkPred.
 *
 * LinkPred: A high performance library for link prediction in complex networks.
 * Copyright (C) 2021 by Said Kerrache.
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
 * @brief Contains the implementation of locally linear graph embedding (LLE).
 */

#ifndef LLE_HPP_
#define LLE_HPP_

#include "LinkPredConfig.hpp"

#ifdef LINKPRED_WITH_ARMADILLO

#include "linkpred/graphalg/encoders/encoder.hpp"

namespace LinkPred {

/**
 * @brief Contains the implementation of locally linear graph embedding (LLE).
 * Reference: Roweis, S. T. and Saul, L. K. (2000). Nonlinear dimensionality
 * reduction by locally linear embedding. Science, 290(5500):2323–2326.
 * @tparam Network The network type.
 */
template<typename Network = UNetwork<>> class LLE: public Encoder<Network> {
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

protected:
	const int DefaultDim = 5; /**< Default embedding dimension. */
	double tol = 0; /**< Tolerance for the eigenvalue problem. */

public:
	/**
	 * Constructor.
	 * @param net The network.
	 */
	LLE(std::shared_ptr<Network const> net) :
			Encoder<Network>(net) {
		dim = DefaultDim;
		name = "LLE";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	LLE(LLE const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	LLE& operator =(LLE const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	LLE(LLE &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	LLE& operator =(LLE &&that) = default;

	/**
	 * Initialize encoder.
	 */
	virtual void init();

	/**
	 * Encode the network.
	 */
	virtual void encode();

	/**
	 * @return The default embedding dimension.
	 */
	const int getDefaultDim() const {
		return DefaultDim;
	}

	/**
	 * @return The tolerance for the eigenvalue problem.
	 */
	double getTol() const {
		return tol;
	}

	/**
	 * Set the tolerance for the eigenvalue problem.
	 * @param tol The tolerance for the eigenvalue problem.
	 */
	void setTol(double tol) {
		this->tol = tol;
	}

	/**
	 * Set the edge weight map.
	 * @param weightMap The edge weight map.
	 */
	virtual void setWeightMap(const WeightMapSP &weightMap) {
		throw std::runtime_error(
				"LLE encoder does not support weighted networks.");
	}

	/**
	 * Destructor.
	 */
	virtual ~LLE() = default;
};

} /* namespace LinkPred */

#endif /* LINKPRED_WITH_ARMADILLO */
#endif /* LLE_HPP_ */
