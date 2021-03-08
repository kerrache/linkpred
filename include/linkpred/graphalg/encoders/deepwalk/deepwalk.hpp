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
 * @brief DeepWalk encoder.
 */

/*
 * Original code license: MIT (attached).
 */

#ifndef DEEPWALK_HPP_
#define DEEPWALK_HPP_

#include "linkpred/graphalg/encoders/encoder.hpp"
#include "linkpred/utils/randomgen.hpp"
#include "linkpred/numerical/mf/fastsig.hpp"
#include <cmath>
#include <stdexcept>
#include <cstddef>

namespace LinkPred {

/**
 * @brief DeepWalk encoder.
 * Reference: Perozzi, B., Al-Rfou, R., and Skiena, S. (2014). Deepwalk: Online learning of social representations.
 * In Proceedings of the 20th ACM SIGKDD International Conference on Knowledge Discovery and Data
 * Mining, KDD ’14, pages 701–710, New York, NY, USA. Association for Computing Machinery.
 * This implementation is based on the code https://github.com/xgfs/deepwalk-c
 * @tparam Network The network type.
 */
template<typename Network = UNetwork<>> class DeepWalk: public Encoder<Network> {
public:
	using typename Encoder<Network>::NodeID; /**< Node ID type. */
	using typename Encoder<Network>::CodeMap; /**< Code map type. */
	using typename Encoder<Network>::CodeMapSP; /**< Shared pointer to a node code map. */
	using typename Encoder<Network>::WeightMap; /**< Edge weight map type. */
	using typename Encoder<Network>::WeightMapSP; /**< Shared pointer to an edge weight map. */

protected:
	using Encoder<Network>::net; /**< The network. */
	using Encoder<Network>::name; /**< The name of the encoder. */
	using Encoder<Network>::weightMap; /**< Weight map. */
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
	int nbWalks = 10; /**< Number of walks per vertex ("\gamma"). */
	int walkLength = 80; /**< DeepWalk parameter "t" = length of the walk. */
	int windowSize = 10; /**< DeepWalk parameter "w" = window size. */
	long long totalSteps = 0; /**< Total number of steps. */
	long long stepInterval = 10000; /**< Number of steps after which the learning rate is updated. */
	int nbNodeWalks = 100; /**< Implementation parameter, number of walks per node in the PageRank initialization. */
	float alpha = 0.85; /**< Probability of continuing random walk in initialization. */
	const int MaxCodeLength = 64; /**< Maximum code length. */
	const int DefaultDim = 10; /**< Default embedding dimension. */

private:
	RandomGen rng; /**< random number generator. */
	long long nbNodes = 0; /**< Number of nodes. */
	long long nbEdges = 0; /**< Number of edges. */
	std::vector<std::size_t> offsets; /**< The ia array of the network CSR representation. */
	std::vector<std::size_t> edges; /**< The ja array of the network CSR representation. */
	std::vector<std::size_t> degrees; /**< Node degrees. */
	std::vector<int> trainOrder; /**< Cache to keep track of order of nodes. */
	std::vector<float> vecEmb; /**< Vertex embedding, aka DeepWalk's \Phi. */
	std::vector<float> vecCtx; /**< Context embedding. */
	std::vector<uint8_t> hsmCodes; /**< HSM codes for each vertex. */
	std::vector<int> hsmPtrs; /*< HSM pointers for each vertex. */
	std::vector<std::size_t> hsmIndPtrs; /**< HSM offsets for each vertex. */
	std::vector<float> hsmWeights; /**< Weights (probabilities) for constructing HSM tree. */
	FastSig<float> sigmoid; /**< Fast sigmoid. */

	/**
	 * Update the embedding, putting wt gradient in wtCache
	 */
	void update(float *ws, float *wt, float *wtCache, float lr,
			const int label);

	/**
	 * Randomly shuffle and array.
	 */
	void shuffle(std::vector<int> &a, int n) { // shuffles the array a of size n
		for (int i = n - 1; i >= 0; i--) {
			int j = rng.getSInt(0, i);
			int temp = a[j];
			a[j] = a[i];
			a[i] = temp;
		}
	}

	/**
	 * Select a random neighbor. If the node has no neighbor, return -1.
	 */
	long long sampleNeighbor(long long i);

	/**
	 * initializes global arrays of HSM from probs array.
	 */
	void initHsm(std::vector<float> &probs);

	/**
	 * Fills the first argument with random walk counts.
	 */
	void estimatePrRW(std::vector<float> &outputs, int nbSamples, float alpha);

public:
	/**
	 * Constructor.
	 * @param net The network.
	 * @param seed Random number generator seed.
	 */
	DeepWalk(std::shared_ptr<Network const> net, long int seed) :
			Encoder<Network>(net), rng(seed) {
		name = "DPW";
		dim = DefaultDim;
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	DeepWalk(DeepWalk const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	DeepWalk& operator =(DeepWalk const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	DeepWalk(DeepWalk &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	DeepWalk& operator =(DeepWalk &&that) = default;

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
	float getInitLR() const {
		return initLR;
	}

	/**
	 * Set the initial learning rate.
	 * @param initLR Initial learning rate.
	 */
	void setInitLR(float initLR) {
		this->initLR = initLR;
	}

	/**
	 * @return Number of walks per vertex ("\gamma").
	 */
	int getNbWalks() const {
		return nbWalks;
	}

	/**
	 * Set the number of walks per vertex ("\gamma").
	 * @param nbWalks Number of walks per vertex ("\gamma").
	 */
	void setNbWalks(int nbWalks) {
		this->nbWalks = nbWalks;
	}

	/**
	 * @return DeepWalk parameter "t" = length of the walk.
	 */
	int getWalkLength() const {
		return walkLength;
	}

	/**
	 * Set the DeepWalk parameter "t" = length of the walk.
	 * @param walkLength DeepWalk parameter "t" = length of the walk.
	 */
	void setWalkLength(int walkLength) {
		this->walkLength = walkLength;
	}

	/**
	 * @return DeepWalk parameter "w" = window size.
	 */
	int getWindowSize() const {
		return windowSize;
	}

	/**
	 * Set the DeepWalk parameter "w" = window size.
	 * @param windowSize DeepWalk parameter "w" = window size.
	 */
	void setWindowSize(int windowSize) {
		this->windowSize = windowSize;
	}

	/**
	 * Set edge weight map.
	 */
	virtual void setWeightMap(const WeightMapSP &weightMap) {
		throw std::runtime_error(
				"DeepWalk encoder does not support weighted networks.");
	}

	/**
	 * @return The probability of continuing random walk in initialization.
	 */
	float getAlpha() const {
		return alpha;
	}

	/**
	 * Set the probability of continuing random walk in initialization.
	 * @param alpha Probability of continuing random walk in initialization.
	 */
	void setAlpha(float alpha) {
		this->alpha = alpha;
	}

	/**
	 * @return Maximum code length.
	 */
	const int getMaxCodeLength() const {
		return MaxCodeLength;
	}

	/**
	 * @return The number of walks per node in the PageRank initialization.
	 */
	int getNbNodeWalks() const {
		return nbNodeWalks;
	}

	/**
	 * Set the number of walks per node in the PageRank initialization.
	 * @param nbNodeWalks The number of walks per node in the PageRank initialization.
	 */
	void setNbNodeWalks(int nbNodeWalks) {
		this->nbNodeWalks = nbNodeWalks;
	}

	/**
	 * @return The number of steps after which the learning rate is updated.
	 */
	long long getStepInterval() const {
		return stepInterval;
	}

	/**
	 * Set the number of steps after which the learning rate is updated.
	 * @param stepInterval Number of steps after which the learning rate is updated.
	 */
	void setStepInterval(long long stepInterval) {
		this->stepInterval = stepInterval;
	}

	/**
	 * @return The total number of steps.
	 */
	long long getTotalSteps() const {
		return totalSteps;
	}

	/**
	 * @return Default embedding dimension.
	 */
	const int getDefaultDim() const {
		return DefaultDim;
	}

	/**
	 * Destructor.
	 */
	virtual ~DeepWalk() = default;

};

} /* namespace LinkPred */

#endif /* DEEPWALK_HPP_ */
