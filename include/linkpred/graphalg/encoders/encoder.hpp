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
 * @brief Contains the interface of a network encoder.
 */

#ifndef ENCODER_HPP_
#define ENCODER_HPP_

#include "LinkPredConfig.hpp"
#include "linkpred/core/unetwork/unetwork.hpp"
#include "linkpred/numerical/linear/vec.hpp"
#ifdef LINKPRED_WITH_OPENMP
#include <omp.h>
#endif
#ifdef LINKPRED_WITH_MPI
#include <mpi.h>
#endif

namespace LinkPred {

/**
 * @brief The interface of a network encoder.
 * @tparam Network The network type.
 */
template<typename Network = UNetwork<>> class Encoder {
public:
	using NodeID = typename Network::NodeID; /**< Node ID type. */
	using Edge = typename Network::Edge; /**< Edge type. */
	using CodeMap = typename Network::template NodeMap<Vec>; /**< Code map type. */
	using CodeMapSP = typename Network::template NodeMapSP<Vec>; /**< Shared pointer to a code map. */
	using WeightMap = typename Network::template EdgeMap<double>; /**< Edge weight map type. */
	using WeightMapSP = typename Network::template EdgeMapSP<double>; /**< Shared pointer to an edge weight map. */

protected:
	std::shared_ptr<Network const> net; /**< The network. */
	std::string name; /**< The name of the encoder. */
	WeightMapSP weightMap; /**< Edge weights. */
	CodeMapSP codeMap; /**< Code map. */
	int dim = 0; /**< Dimension of the embedding. */
#ifdef LINKPRED_WITH_OPENMP
	bool parallel = false; /**< Enable/disable shared-memory parallelism. */
#endif
#ifdef LINKPRED_WITH_MPI
	bool distributed = false; /**< Enable/disable distributed parallelism. */
	MPI_Comm comm = MPI_COMM_WORLD; /**< The MPI communicator. */
#endif

public:
	/**
	 * Constructor.
	 * @param net The network.
	 */
	Encoder(std::shared_ptr<Network const> net) :
			net(net) {
	}

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	Encoder(Encoder const &that) = default;

	/**
	 * Copy assignment operator.
	 * @param that The object to copy.
	 */
	Encoder& operator =(Encoder const &that) = default;

	/**
	 * Move constructor.
	 * @param that The object to move.
	 */
	Encoder(Encoder &&that) = default;

	/**
	 * Move assignment operator.
	 * @param that The object to move.
	 */
	Encoder& operator =(Encoder &&that) = default;

#ifdef LINKPRED_WITH_OPENMP
	/**
	 * @return Whether shared memory parallelism is enabled.
	 */
	bool isParallel() const {
		return parallel;
	}

	/**
	 * Enable/disable shared memory parallelism.
	 * @param parallel True to enable parallelism, false to disable it.
	 */
	void setParallel(bool parallel) {
		this->parallel = parallel;
	}
#endif

#ifdef LINKPRED_WITH_MPI
	/**
	 * @return Whether distributed memory parallelism is enabled.
	 */
	bool isDistributed() const {
		return distributed;
	}

	/**
	 * Enable/disable distributed memory parallelism.
	 * @param distributed True to enable distributed memory parallelism, false to disable it.
	 */
	void setDistributed(bool distributed) {
		this->distributed = distributed;
	}

	/**
	 * @return The MPI communicator.
	 */
	MPI_Comm getComm() const {
		return comm;
	}

	/**
	 * Set the MPI communicator.
	 * @param comm The new MPI communicator.
	 */
	void setComm(MPI_Comm const & comm) {
		this->comm = comm;
	}
#endif

	/**
	 * Initialize encoder.
	 */
	virtual void init() = 0;

	/**
	 * Encode the network.
	 */
	virtual void encode() = 0;

	/**
	 * Return the encoding of all nodes in the network.
	 * @return The network encoding in the form of a node map.
	 */
	CodeMapSP getNodeCodeMap() {
		return codeMap;
	}

	/**
	 * Return the code of given node.
	 * @param i The node ID.
	 * @return The encoding of node i.
	 */
	Vec getNodeCode(NodeID const &i) {
		return codeMap->at(i);
	}

	/**
	 * Return the code of given edge. This is a default implementation that simply concatenates the two node codes.
	 * Sub-classes can override this implementation and provide a more complex one.
	 * @param i The start node ID.
	 * @param j The end node ID.
	 * @return The encoding of edge (i,j).
	 */
	virtual Vec getEdgeCode(NodeID const &i, NodeID const &j) {
		return Vec(codeMap->at(i), codeMap->at(j));
	}

	/**
	 * Shortcut for the previous method.
	 * @param e The edge.
	 * @return The encoding of e.
	 */
	Vec getEdgeCode(Edge const &e) {
		return this->getEdgeCode(net->start(e), net->end(e));
	}

	/**
	 * This corresponds to the default edge encoding (see the method getEdgeCode).
	 * @return The dimension of the edge codes.
	 */
	virtual int getEdgeCodeDim() const {
		return 2 * dim;
	}

	/**
	 * @return The dimension of the embedding (node embedding).
	 */
	int getDim() const {
		return dim;
	}

	/**
	 * Set the dimension of the embedding (node embedding).
	 * @param dim The dimension of the embedding (node embedding).
	 */
	void setDim(int dim) {
		this->dim = dim;
	}

	/**
	 * @return The name of the encoder.
	 */
	const std::string& getName() const {
		return name;
	}

	/**
	 * Set the name of the encoder.
	 * @param name The new name of the encoder.
	 */
	void setName(const std::string &name) {
		this->name = name;
	}

	/**
	 * @return The edge weight map.
	 */
	const WeightMapSP& getWeightMap() const {
		return weightMap;
	}

	/**
	 * Set The edge weight map.
	 * @param weightMap The edge weight map.
	 */
	virtual void setWeightMap(const WeightMapSP &weightMap) {
		this->weightMap = weightMap;
	}

	/**
	 * Destructor.
	 */
	virtual ~Encoder() = default;

};

} /* namespace LinkPred */

#endif /* ENCODER_HPP_ */
