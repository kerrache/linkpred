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
 * @brief The LINE (Large-scale Information Network Embedding) encoder.
 */

/*
 * Original code license: Apache 2.0 (attached).
 */

#ifndef LINE_HPP_
#define LINE_HPP_

#include "linkpred/graphalg/encoders/encoder.hpp"
#include "linkpred/utils/randomgen.hpp"
#include "linkpred/numerical/mf/fastsig.hpp"
#include <cmath>
#include <stdexcept>

namespace LinkPred {

/**
 * @brief LINE encoder.
 *
 * This implementation is based on the code https://github.com/tangjianpku/LINE
 * Reference: Tang, J., Qu, M., Wang, M., Zhang, M., Yan, J., and Mei, Q. (2015).
 * LINE: Large-Scale Information Network Embedding, page 1067–1077. International
 * World Wide Web Conferences Steering Committee, Republic and Canton of Geneva, CHE.
 * @tparam Network The network type.
 */
template<typename Network = UNetwork<>> class LINE: public Encoder<Network> {
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

	float initLR = 0.025f; /**< Initial learning rate. */
	int nbNegSamples = 5; /**< Number of negative samples. */
	int order = 1; /**< Order of similarity. Possible values are: 1 (first order), 2 (second order), and 12 (concatenate first and second order). When order = 12, the dimension is split between first order codes and second order codes. */
	long long negSize = 1e8; /**< Size of the table used in negative sampling. */
	float negSamplingPower = 0.75; /**< Power used in negative sampling. */
	long long stepInterval = 10000; /**< Number of steps after which the learning rate is updated. */
	int maxDepth = 2; /**< The maximum depth in the Breadth-First-Search. */
	int threshold = 10; /**< For vertex whose degree is less than threshold, we will expand its neighbors until the degree reaches.*/
	const int DefaultDim = 10; /**< Default embedding dimension. */
	long long totalSteps = 0; /**< Total number of steps. */
	bool enableReconstruct = false; /**< Whether to expand the network. */

private:
	RandomGen rng; /**< random number generator. */
	long long nbNodes = 0; /**< Number of nodes. */
	long long nbEdges = 0; /**< Number of edges. */
	std::vector<int> edgeFrom;
	std::vector<int> edgeTo;
	std::vector<double> wDegree; /**< Weighted degree (equals the degree in an unweighted network). */
	std::vector<double> edgeWeight;
	std::vector<long long> alias;
	std::vector<double> prob;
	std::vector<int> negTable;
	std::vector<float> vecEmb; /**< Vertex embedding, aka DeepWalk's \Phi. */
	std::vector<float> vecCtx; /**< Context embedding. */
	FastSig<float> sigmoid; /**< Fast sigmoid. */
	std::shared_ptr<Network> recNet; /**< The reconstructed network. */
	WeightMapSP recWeightMap; /**< Edge weights in the reconstructed network. */

	/**
	 * A helper function which allows to initialize data either from the original or
	 * the reconstructed network.
	 * @param initNet The network used for initialization.
	 * @param initWeightMap The edge weight map used for initialization.
	 */
	void init(std::shared_ptr<Network const> initNet,
			WeightMapSP initWeightMap);

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

	/**
	 * Gradient descent update.
	 */
	void update(float *vecU, float *vecV, float *vecErr, float lr,
			const int label);

	/**
	 * The reconstruction step. This is more of a misnomer as the network is actually "extended" by adding more
	 * edges; No reconstruction actually takes place.
	 */
	void reconstruct();

	/**
	 * Encode with a given order.
	 */
	void encode(int encOrder);

public:

	/**
	 * Constructor.
	 * @param net The network.
	 * @param seed Random number generator seed.
	 */
	LINE(std::shared_ptr<Network const> net, long int seed) :
			Encoder<Network>(net), rng(seed) {
		name = "LIN";
		dim = DefaultDim;
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	LINE(LINE const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	LINE& operator =(LINE const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	LINE(LINE &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	LINE& operator =(LINE &&that) = default;

	/**
	 * Initialize encoder.
	 */
	virtual void init();

	/**
	 * Encode the network.
	 */
	virtual void encode();

	/**
	 * @return Initial learning rate.
	 */
	float getInitLr() const {
		return initLR;
	}

	/**
	 * Set the initial learning rate.
	 * @param initLR The initial learning rate.
	 */
	void setInitLr(float initLR) {
		this->initLR = initLR;
	}

	/**
	 * @return Number of negative samples.
	 */
	int getNbNegSamples() const {
		return nbNegSamples;
	}

	/**
	 * Set the number of negative samples.
	 * @param nbNegSamples The number of negative samples.
	 */
	void setNbNegSamples(int nbNegSamples) {
		this->nbNegSamples = nbNegSamples;
	}

	/**
	 * @return The order of similarity. Possible values are: 1 (first order), 2 (second order),
	 * and 12 (concatenate first and second order). When order = 12, the embedding dimension
	 * is split between first order codes and second order codes.
	 */
	int getOrder() const {
		return order;
	}

	/**
	 * Set the order of similarity. Possible values are: 1 (first order), 2 (second order),
	 * and 12 (concatenate first and second order). When order = 12, the embedding dimension
	 * is split between first order codes and second order codes.
	 * @param order Order of similarity.
	 */
	void setOrder(int order) {
		if (order != 1 && order != 2 && order != 12) {
			throw std::runtime_error(
					"The parameter order can only take the values 1, 2, and 12");
		}
		this->order = order;
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
	 * @return Whether the network is reconstructed (expanded).
	 */
	bool isEnableReconstruct() const {
		return enableReconstruct;
	}

	/**
	 * Set whether the network is reconstructed (expanded).
	 * @param enableReconstruct Whether the network is reconstructed (expanded).
	 */
	void setEnableReconstruct(bool enableReconstruct) {
		this->enableReconstruct = enableReconstruct;
	}

	/**
	 * @return The reconstructed network.
	 */
	const std::shared_ptr<const Network> getRecNet() const {
		return recNet;
	}

	/**
	 * @return The weight map of the reconstructed network.
	 */
	const WeightMapSP& getRecWeightMap() const {
		return recWeightMap;
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
	virtual ~LINE() = default;
};

} /* namespace LinkPred */

#endif /* LINE_HPP_ */
