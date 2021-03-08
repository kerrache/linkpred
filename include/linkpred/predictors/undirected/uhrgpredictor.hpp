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

// Original copyright notice
/*
 * fitHRG - fits a hierarchical random graph (hrg) model to data
 * Copyright (C) 2005-2012 Aaron Clauset
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * See http: *www.gnu.org/licenses/gpl.txt for more details.
 *
 * ****************************************************************************************************
 * Author       : Aaron Clauset  ( aaronc@santafe.edu | http: *www.santafe.edu/~aaronc/ )
 * Collaborators: Cristopher Moore and Mark Newman
 * Project      : Hierarchical Random Graphs
 * Location     : University of New Mexico, Dept. of Computer Science AND Santa Fe Institute
 * Created      : 21 June 2006
 * Modified     : 16 June 2007
 *			 : 26 December 2007 (cleaned up for public consumption)
 *			 : 27 May 2012      (workaround for quicksort bug)
 *
 * ****************************************************************************************************
 *
 * This program runs the MCMC with HRG model and the input graph G, samples the likelihoods of edges
 * that are not present in G, and writes out a list of them rank-ordered by their average likelihood
 * under the sampled HRG models. The program includes a convergence criterion that attempts to detect
 * when the MCMC has converged on the equilibrium set.
 *
 * ****************************************************************************************************
 *
 *  See http: *www.santafe.edu/~aaronc/randomgraphs/ for more information, updates to the code, etc.
 */

/**
 * \file
 * @ingroup Predictors
 * @brief Contains the implementation of the HRG link predictor.
 */

#ifndef UHRGPREDICTOR_HPP_
#define UHRGPREDICTOR_HPP_

#include <linkpred/predictors/undirected/uhrgpredictor/dendro_pr.hpp>
#include <linkpred/predictors/undirected/uhrgpredictor/graph_pr.h>
#include <linkpred/predictors/undirected/uhrgpredictor/graph_simp.h>
#include <linkpred/predictors/undirected/uhrgpredictor/rbtree.h>
#include <linkpred/predictors/undirected/ulpredictor.hpp>
#include "linkpred/utils/randomgen.hpp"
#include <string>
#include <map>

namespace LinkPred {

/**
 * @brief HRG predictor.
 * @details This is actually a modified and wrapped version the code provided by the authors.
 * @tparam Network The network type.
 * @tparam EdgeRndIt A random iterator type used to iterate on edges.
 * @tparam ScoreRndIt A random iterator type used to iterate on scores.
 */
template<typename Network = UNetwork<>,
		typename EdgeRndIt = typename std::vector<typename Network::Edge>::const_iterator,
		typename ScoreRndIt = typename std::vector<double>::iterator,
		typename EdgeRndOutIt = typename std::vector<typename Network::Edge>::iterator> class UHRGPredictor: public ULPredictor<
		Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt> {

	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::net; /**< The network. */
#ifdef LINKPRED_WITH_OPENMP
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::parallel; /**< Whether the predictor runs in parallel. */
#endif
	using ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::name; /**< The name of the predictor. */
	using NodeID = typename ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::NodeID; /**< The node IDs type. */
	using Edge = typename ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>::Edge; /**< The edges type. */

	/**
	 * @brief A p-block structure.
	 */
	struct pblock {
		double L;
		int i;
		int j;
	};

	/**
	 * @brief I/O parameters.
	 */
	struct ioparameters {
		int n; /**< Number vertices in input graph. */
		int m; /**< Number of edges in input graph. */
		std::string d_dir; /**< Working directory. */
		std::string f_pairs; /**< (in-file) original adjacency list of graph (from researcher). */
		std::string f_hrg; /**< (in-file) dendrogram seed file for MCMC (from researcher). */
		std::string f_X1; /**< (out-file) list of non-edges, ranked by model-averaged likelihood. */
		std::string s_scratch; /**< Scratch space for building filenames. */
		std::string s_tag; /**< User defined filename tag. */
		int timer; /**< Timer for reading input. */
		bool flag_timer; /**< (flag) for timer. */
		bool flag_compact; /**< (flag) T: compress number of trials. */
		bool flag_f_hrg; /**< (flag) T: use f_hrg file as seed. */
		std::string start_time; /**< time simulation was started. */
	};

protected:
	ioparameters ioparm; /**< Program parameters. */
	dendro* d = nullptr; /**< Inferred dendrograph data structure. */
	simpleGraph* g = nullptr; /**< Base graph read from file. */
	int nbSamples = 10000; /**< Number of samples to take for predictions. */
	int nbBeans = 25; /**< Number of bins in edge statistics histogram. */
	long int seed; /**< The random number generator's seed. */
	RandomGen rng; /**< Mersenne Twister random number generator instance. */
	std::map<Edge, double> scores; /**< Links scores. */

	/**
	 * Run the MCMC until we've found equilibrium.
	 * @return True if no error occurs.
	 */
	bool mcmcEquilibriumFind();

	/**
	 * Sample likelihoods for missing connections.
	 * @return True if no error occurs.
	 */
	bool mcmcEquilibriumSample();

	/**
	 * Convert the network into the format used by HRG.
	 * @return True if conversion is successful.
	 */
	bool toHRGNet();

	/**
	 * Setup the predictor.
	 */
	void setup();

	/**
	 * Cleanup at the end.
	 */
	void cleanup();

public:

	/**
	 * @param net The network.
	 * @param seed The random number generator's seed.
	 */
	UHRGPredictor(std::shared_ptr<Network const> net, long int seed) :
			ULPredictor<Network, EdgeRndIt, ScoreRndIt, EdgeRndOutIt>(net), seed(
					seed), rng(seed) {
		name = "HRG";
	}

	/**
	 * Copy constructor.
	 */
	UHRGPredictor(UHRGPredictor const & that) = delete;

	/**
	 * Copy assignment operator.
	 */
	UHRGPredictor & operator =(UHRGPredictor const & that) = delete;

	/**
	 * Move constructor.
	 */
	UHRGPredictor(UHRGPredictor && that) = delete;

	/**
	 * Move assignment operator.
	 */
	UHRGPredictor & operator =(UHRGPredictor && that) = delete;

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
		return scores.at(e);
	}

	/**
	 * Predict the links.
	 * @param begin Beginning of the links to be predicted.
	 * @param end end of the links to be predicted.
	 * @param scores Beginning of scores.
	 */
	virtual void predict(EdgeRndIt begin, EdgeRndIt end, ScoreRndIt scores);

	/**
	 * @return The number of bins.
	 */
	int getNbBeans() const {
		return nbBeans;
	}

	/**
	 * Set the number of bins.
	 * @param nbBeans The new number of bins.
	 */
	void setNbBeans(int nbBeans) {
		this->nbBeans = nbBeans;
	}

	/**
	 * @return The number of samples.
	 */
	int getNbSamples() const {
		return nbSamples;
	}

	/**
	 * Set the number of samples.
	 * @param nbSamples The new number of samples.
	 */
	void setNbSamples(int nbSamples) {
		this->nbSamples = nbSamples;
	}

	/**
	 * Destructor.
	 */
	virtual ~UHRGPredictor();

};

} /* namespace LinkPred */

#endif /* UHRGPREDICTOR_HPP_ */
