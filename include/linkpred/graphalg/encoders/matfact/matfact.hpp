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
 * @brief Contains the implementation of an algorithm for embedding a network
 * using matrix factorization.
 * Reference:
 * Koren, Y., Bell, R., and Volinsky, C. (2009). Matrix factorization techniques
 * for recommender systems. Computer, 42(8):30–37
 * Ahmed, A., Shervashidze, N., Narayanamurthy, S., Josifovski, V., and Smola,
 * A. J. (2013). Distributed large-scale natural graph factorization.
 * In Proceedings of the 22nd International Conference on World Wide Web, WWW ’13,
 * pages 37–48, New York, NY, USA. Association for Computing Machinery.
 */

#ifndef MATFACT_HPP_
#define MATFACT_HPP_

#include "linkpred/graphalg/encoders/encoder.hpp"
#include "linkpred/graphalg/encoders/matfact/matfactcg.hpp"
#include "linkpred/utils/randomgen.hpp"
#include <cmath>

namespace LinkPred {

/**
 * @brief Contains the implementation of an algorithm for embedding a network
 * using matrix factorization.
 * Reference:
 * Koren, Y., Bell, R., and Volinsky, C. (2009). Matrix factorization techniques
 * for recommender systems. Computer, 42(8):30–37
 * Ahmed, A., Shervashidze, N., Narayanamurthy, S., Josifovski, V., and Smola,
 * A. J. (2013). Distributed large-scale natural graph factorization.
 * In Proceedings of the 22nd International Conference on World Wide Web, WWW ’13,
 * pages 37–48, New York, NY, USA. Association for Computing Machinery.
 * @tparam Network The network type.
 */
template<typename Network = UNetwork<>> class MatFact: public Encoder<Network> {
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
	double lambda = 0.001; /**< Regularization coeffcient. */
	double posRatio = 1.0; /**< Ratio of positive edges used for fitting the model. */
	double negRatio = 0.5; /**< Ratio of negative edges used for fitting the model. */
	double tol = 1.0e-8; /**< Tolerance for stopping the optimization. */

private:
	RandomGen rng; /**< Random number generator. */

public:
	/**
	 * Constructor.
	 * @param net The network.
	 * @param seed Random number generator seed.
	 */
	MatFact(std::shared_ptr<Network const> net, long int seed) :
			Encoder<Network>(net), rng(seed) {
		dim = DefaultDim;
		name = "MFC";
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	MatFact(MatFact const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	MatFact& operator =(MatFact const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	MatFact(MatFact &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	MatFact& operator =(MatFact &&that) = default;

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
	 * @return The regularization coeffcient.
	 */
	double getLambda() const {
		return lambda;
	}

	/**
	 * Set the regularization coeffcient.
	 * @param lambda The regularization coeffcient.
	 */
	void setLambda(double lambda) {
		this->lambda = lambda;
	}

	/**
	 * @return The ratio of negative edges used for fitting the model.
	 */
	double getNegRatio() const {
		return negRatio;
	}

	/**
	 * Set the ratio of negative edges used for fitting the model.
	 * @param negRatio The ratio of negative edges used for fitting the model.
	 */
	void setNegRatio(double negRatio) {
		this->negRatio = negRatio;
	}

	/**
	 * @return The ratio of positive edges used for fitting the model.
	 */
	double getPosRatio() const {
		return posRatio;
	}

	/**
	 * Set the ratio of positive edges used for fitting the model.
	 * @param posRatio The ratio of positive edges used for fitting the model.
	 */
	void setPosRatio(double posRatio) {
		this->posRatio = posRatio;
	}

	/**
	 * @return The tolerance for stopping the optimization.
	 */
	double getTol() const {
		return tol;
	}

	/**
	 * Set the tolerance for stopping the optimization.
	 * @param tol The tolerance for stopping the optimization.
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
				"MatFact encoder does not support weighted networks.");
	}

	/**
	 * Destructor.
	 */
	virtual ~MatFact() = default;
};

} /* namespace LinkPred */

#endif /* MATFACT_HPP_ */
