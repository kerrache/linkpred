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
 * @brief A node2vec encoder.
 */

/*
 * Original code license: MIT (attached).
 */

#ifndef NODE2VEC_HPP_
#define NODE2VEC_HPP_

#include "linkpred/graphalg/encoders/encoder.hpp"
#include "linkpred/utils/randomgen.hpp"
#include "linkpred/numerical/mf/fastsig.hpp"
#include <cmath>
#include <stdexcept>

namespace LinkPred {

/**
 * @brief Node2Vec encoder.
 * References: Grover, A. and Leskovec, J. (2016). Node2vec: Scalable feature learning
 * for networks. In Proceedings of the 22nd ACM SIGKDD International Conference on
 * Knowledge Discovery and Data Mining, KDD’16, pages 855–864, New York, NY, USA.
 * Association for Computing Machinery.
 * This implementation is based on the code https://github.com/xgfs/node2vec-c
 * @tparam Network The network type.
 */
template<typename Network = UNetwork<>> class Node2Vec: public Encoder<Network> {
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
	int nbWalks = 10; /**< DeepWalk parameter "\gamma" = walks per vertex. */
	int walkLength = 80; /**< DeepWalk parameter "t" = length of the walk. */
	int windowSize = 10; /**< DeepWalk parameter "w" = window size. */
	int nbNegSamples = 5; /**< Number of negative samples. */
	float p = 1; /**< The parameter p. */
	float q = 1; /**< The parameter q. */
	long long totalSteps = 0; /**< Total number of steps. */
	long long stepInterval = 10000; /**< Number of steps after which the learning rate is updated. */
	float subSample = 0.0; /**< Sub-sample size. */
	const int DefaultDim = 10; /**< Default embedding dimension. */

private:
	RandomGen rng; /**< random number generator. */
	long long nbNodes = 0; /**< Number of nodes. */
	long long nbEdges = 0; /**< Number of edges. */
	std::vector<std::size_t> offsets; /**< The ia array of the network CSR representation. */
	std::vector<std::size_t> edges; /**< The ja array of the network CSR representation. */
	std::vector<std::size_t> degrees; /**< Node degrees. */
	std::vector<int> trainOrder; /**< Cache to keep track of order of nodes. */
	std::vector<unsigned long long> edgeOffsets; /**< Alias table pointers. */
	std::vector<int> n2vJs;
	std::vector<float> n2vQs;
	std::vector<int> negJs;
	std::vector<float> negQs;
	std::vector<float> nodeCnts;
	std::vector<float> vecEmb; /**< Vertex embedding, aka DeepWalk's \Phi. */
	std::vector<float> vecCtx; /**< Context embedding. */
	FastSig<float> sigmoid; /**< Fast sigmoid. */

	/**
	 * Update the embedding, putting wt gradient in wtCache
	 */
	void update(float *ws, float *wt, float *wtCache, float lr,
			const int label);

	/**
	 * Initialize walker. Assumes probs are normalized.
	 */
	void initWalker(int n, int *j, float *probs);

	/**
	 * Random selection.
	 */
	int walkerDraw(const int n, float *q, int *j) {
		int kk = rng.getSInt(0, n - 1);
		return rng.getDouble(0, 1) < q[kk] ? kk : j[kk];
	}

	/**
	 * Check if there is an edge.
	 */
	bool hasEdge(int from, int to) {
		return std::binary_search(&edges[offsets[from]],
				&edges[offsets[from + 1]], to);
	}

	/**
	 * shuffles the array a of size n.
	 */
	void shuffle(std::vector<int> &a, int n) {
		for (int i = n - 1; i >= 0; i--) {
			int j = rng.getSInt(0, i);
			int temp = a[j];
			a[j] = a[i];
			a[i] = temp;
		}
	}

public:
	/**
	 * Constructor.
	 * @param net The network.
	 * @param seed Random number generator seed.
	 */
	Node2Vec(std::shared_ptr<Network const> net, long int seed) :
			Encoder<Network>(net), rng(seed) {
		name = "N2V";
		dim = DefaultDim;
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	Node2Vec(Node2Vec const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	Node2Vec& operator =(Node2Vec const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	Node2Vec(Node2Vec &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	Node2Vec& operator =(Node2Vec &&that) = default;

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
	 * @return The parameter p.
	 */
	float getP() const {
		return p;
	}

	/**
	 * Set the parameter p.
	 * @param p The parameter p.
	 */
	void setP(float p) {
		this->p = p;
	}

	/**
	 * @return The parameter q.
	 */
	float getQ() const {
		return q;
	}

	/**
	 * Set the parameter q.
	 * @param q The parameter q.
	 */
	void setQ(float q) {
		this->q = q;
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
	virtual void setWeightMap(const WeightMapSP& weightMap) {
		throw std::runtime_error(
				"Node2Vec encoder does not support weighted networks.");
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
	 * @return Sub-sample size.
	 */
	float getSubSample() const {
		return subSample;
	}

	/**
	 * Set sub-sample size.
	 * @param subSample Sub-sample size.
	 */
	void setSubSample(float subSample) {
		this->subSample = subSample;
	}

	/**
	 * Destructor.
	 */
	virtual ~Node2Vec() = default;

};

} /* namespace LinkPred */

#endif /* NODE2VEC_HPP_ */
