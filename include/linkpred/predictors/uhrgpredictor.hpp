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
 * @brief Contains the implementation of the HRG link predictor.
 */

#ifndef UHRGPREDICTOR_HPP_
#define UHRGPREDICTOR_HPP_

#include <linkpred/predictors/uhrgpredictor/dendro_pr.hpp>
#include <linkpred/predictors/uhrgpredictor/graph_pr.h>
#include <linkpred/predictors/uhrgpredictor/graph_simp.h>
#include <linkpred/predictors/uhrgpredictor/rbtree.h>
#include <linkpred/predictors/ulpredictor.hpp>
#include "linkpred/utils/randomgen.hpp"
#include <string>
#include <map>

namespace LinkPred {

/**
 * @brief HRG predictor.
 * @details This is actually a modified and wrapped version the code provided by the authors.
 * @tparam NetworkT The network type.
 * @tparam EdgesRandomIteratorT A random iterator type used to iterate on edges.
 * @tparam ScoresRandomIteratorT A random iterator type used to iterate on scores.
 */
template<typename NetworkT = UNetwork<>,
		typename EdgesRandomIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::const_iterator,
		typename ScoresRandomIteratorT = typename std::vector<double>::iterator,
		typename EdgesRandomOutputIteratorT = typename std::vector<
				typename NetworkT::EdgeType>::iterator> class UHRGPredictor: public ULPredictor<
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
	int numSamples = 10000; /**< Number of samples to take for predictions. */
	int numBins = 25; /**< Number of bins in edge statistics histogram. */
	long int seed; /**< The random number generator's seed. */
	RandomGen rng; /**< Mersenne Twister random number generator instance. */
	std::map<EdgeType, double> scores; /**< Links scores. */

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
	UHRGPredictor(std::shared_ptr<NetworkT const> net, long int seed) :
			ULPredictor<NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
					EdgesRandomOutputIteratorT>(net), seed(seed), rng(seed) {
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
	virtual double score(EdgeType const & e) {
		return scores.at(e);
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
	 * @return The number of bins.
	 */
	int getNumBins() const {
		return numBins;
	}

	/**
	 * Set the number of bins.
	 * @param numBins The new number of bins.
	 */
	void setNumBins(int numBins) {
		this->numBins = numBins;
	}

	/**
	 * @return The number of samples.
	 */
	int getNumSamples() const {
		return numSamples;
	}

	/**
	 * Set the number of samples.
	 * @param numSamples The new number of samples.
	 */
	void setNumSamples(int numSamples) {
		this->numSamples = numSamples;
	}

	/**
	 * Destructor.
	 */
	virtual ~UHRGPredictor();

};

}
/* namespace LinkPred */

#endif /* UHRGPREDICTOR_HPP_ */
