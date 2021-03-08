/*
 * This file is part of LinkPred.
 *
 * LinkPred: A high performance library for link prediction in complex networks.
 * Copyright (C) 2021  by Said Kerrache.
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
 * @brief LargeVis encoder.
 */

/*
 * Original code license: Apache 2.0 (attached).
 */

#ifndef LARGEVIS_HPP_
#define LARGEVIS_HPP_

#include "linkpred/graphalg/encoders/encoder.hpp"
#include "linkpred/utils/randomgen.hpp"
#include <string>
#include <vector>

namespace LinkPred {

/**
 * @brief LargeVis encoder.
 * Reference: Tang, J., Liu, J., Zhang, M., and Mei, Q. (2016b). Visualizing large-scale
 * and high-dimensional data. In Bourdeau, J., Hendler, J., Nkambou, R., Horrocks, I.,
 * and Zhao, B. Y., editors, WWW, pages 287–297. ACM.
 * This implementation is based on the code https://github.com/lferry007/LargeVis
 * @tparam Network The network type.
 */
template<typename Network = UNetwork<>> class LargeVis: public Encoder<Network> {
public:
	using typename Encoder<Network>::NodeID; /**< Node ID type. */
	using typename Encoder<Network>::CodeMap; /**< Code map type. */
	using typename Encoder<Network>::CodeMapSP; /**< Shared pointer to a node code map. */
	using typename Encoder<Network>::WeightMap; /**< Edge weight map type. */
	using typename Encoder<Network>::WeightMapSP; /**< Shared pointer to an edge weight map. */

protected:
	using Encoder<Network>::net; /**< The network. */
	using Encoder<Network>::name; /**< The name of the encoder. */
	using Encoder<Network>::weightMap; /**< Code map. */
	using Encoder<Network>::codeMap; /**< Code map. */
	using Encoder<Network>::dim; /**< Dimension of the embedding. */
#ifdef LINKPRED_WITH_OPENMP
	using Encoder<Network>::parallel; /**< Enable/disable shared-memory parallelism. */
#endif
#ifdef LINKPRED_WITH_MPI
	using Encoder<Network>::distributed; /**< Enable/disable distributed parallelism. */
	using Encoder<Network>::comm; /**< The MPI communicator. */
#endif
	long long nbSamples = 10000; /**< Nnumber of samples. */
	long long nbNegSamples = 5; /**< Number of negative samples. */
	float initLR = 1; /**< Initial learning rate. */
	float gamma = 7; /**< The parameter gamma. */
	long long negSize = 1e8; /**< Size of the table used in negative sampling. */
	float negSamplingPower = 0.75; /**< Power used in negative sampling. */
	const int DefaultDim = 10; /**< Default embedding dimension. */

private:
	RandomGen rng; /**< random number generator. */
	long long nbNodes = 0;
	long long nbEdges = 0;
	std::vector<float> vecEmb;
	std::vector<long long> head;
	std::vector<long long> next;
	std::vector<long long> reverse;
	std::vector<int> edgeFrom;
	std::vector<int> edgeTo;
	std::vector<float> edgeWeight;
	std::vector<int> negTable;
	std::vector<long long> alias;
	std::vector<float> prob;

	/**
	 * Initialize alias table.
	 */
	void initAliasTable();

	/**
	 * Sample an edge randomly.
	 */
	long long sampleEdge();

	/**
	 * Initialize negTable.
	 */
	void initNegTable();

public:
	/**
	 * Constructor.
	 * @param net The network.
	 * @param seed Random number generator seed.
	 */
	LargeVis(std::shared_ptr<Network const> net, long int seed) :
			Encoder<Network>(net), rng(seed) {
		name = "LVS";
		dim = DefaultDim;
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	LargeVis(LargeVis const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	LargeVis& operator =(LargeVis const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	LargeVis(LargeVis &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	LargeVis& operator =(LargeVis &&that) = default;

	/**
	 * Initialize encoder.
	 */
	virtual void init();

	/**
	 *
	 * Encode the network.
	 */
	virtual void encode();

	/**
	 * Set nbSamples using a simple rule of thumb.
	 */
	void setNbSamples();

	/**
	 * @return The parameter gamma.
	 */
	float getGamma() const {
		return gamma;
	}

	/**
	 * Set the parameter gamma.
	 * @param gamma The parameter gamma.
	 */
	void setGamma(float gamma) {
		this->gamma = gamma;
	}

	/**
	 * @return Initial learning rate.
	 */
	float getInitLR() const {
		return initLR;
	}

	/**
	 * Set the initial learning rate.
	 * @param initLR The initial learning rate.
	 */
	void setInitLR(float initLR) {
		this->initLR = initLR;
	}

	/**
	 * @return Number of negative samples.
	 */
	long long getNbNegSamples() const {
		return nbNegSamples;
	}

	/**
	 * Set the number of negative samples.
	 * @param nbNegSamples The number of negative samples.
	 */
	void setNbNegSamples(long long nbNegSamples) {
		this->nbNegSamples = nbNegSamples;
	}

	/**
	 * @return Number of samples.
	 */
	long long getNbSamples() const {
		return nbSamples;
	}

	/**
	 * Set the number of samples.
	 * @param nbSamples The number of samples.
	 */
	void setNbSamples(long long nbSamples) {
		this->nbSamples = nbSamples;
	}

	/**
	 * @return The size of the table used in negative sampling.
	 */
	long long getNegSize() const {
		return negSize;
	}

	/**
	 * Set the size of the table used in negative sampling.
	 * @param negSize Size of the table used in negative sampling.
	 */
	void setNegSize(long long negSize) {
		this->negSize = negSize;
	}

	/**
	 * @param Power used in negative sampling.
	 */
	float getNegSamplingPower() const {
		return negSamplingPower;
	}

	/**
	 * Set the power used in negative sampling.
	 * @param negSamplingPower Power used in negative sampling.
	 */
	void setNegSamplingPower(float negSamplingPower) {
		this->negSamplingPower = negSamplingPower;
	}

	/**
	 * @return The default embedding dimension.
	 */
	const int getDefaultDim() const {
		return DefaultDim;
	}

	/**
	 * Destructor.
	 */
	virtual ~LargeVis() = default;

};

} /* namespace LinkPred */

#endif /* LARGEVIS_HPP_ */
