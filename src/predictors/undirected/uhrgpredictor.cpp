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

#include <linkpred/predictors/undirected/uhrgpredictor.hpp>
#include "linkpred/utils/log.hpp"

namespace LinkPred {


template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UHRGPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::init() {
	cleanup();
	ioparm.timer = 20;				//
	ioparm.s_tag = "";				//
	ioparm.f_hrg = "";				//
	ioparm.flag_f_hrg = false;			//
	scores.clear();
	d = new dendro(rng.getInt()); // create hrg data structure
	toHRGNet();
	d->buildDendrogram(); // setup the data structure for the test
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UHRGPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::learn() {

	logger(logDebug1, "Beginning convergence to equilibrium...")
	if (!(mcmcEquilibriumFind())) {
		throw std::runtime_error("No MCMC convergence.");
	}	// run it to equilibrium
	logger(logDebug1, "Done")

	logger(logDebug1, "Beginning sampling...")
	if (!(mcmcEquilibriumSample())) {
		throw std::runtime_error("Error when sampling.");
	}	// sample likelihoods for missing connections
	logger(logDebug1, "Done")

	for (auto it = net->nonEdgesBegin(); it != net->nonEdgesEnd(); ++it) {
		auto i = Network::start(*it);
		auto j = Network::end(*it);
		double w = g->getAdjacency(i, j);
		if ((w != 0) && (w != 1)) {
			throw std::logic_error("The adjacency matrix can only have 0 or 1");
		}
		if (w == 1) {
			throw std::logic_error(
					"The adjacency matrix entry for a non-existing link must be 0");
		}
		double temp = d->g->getAdjacencyAverage(i, j);
		scores[*it] = temp * (1.0 + rng.getDouble(0, 1) / 1000.0);
	}
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UHRGPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::predict(EdgeRndIt begin,
		EdgeRndIt end, ScoreRndIt oscores) {
	logger(logDebug, "Predicting links...")
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		*(oscores + (it - begin)) = scores.at(*it);
	}
	logger(logDebug, "Done")
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> UHRGPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::~UHRGPredictor() {
	cleanup();
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> bool UHRGPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::mcmcEquilibriumFind() {
	double dL, likeli, bestL, oldMeanL, newMeanL;
	bool flagTaken, flagEq;
	int t = 1;

	// We want to run the MCMC until we've found equilibrium; we use the heuristic of the
	// average log-likelihood (which is exactly the entropy) over X steps being very close
	// to the average log-likelihood (entropy) over the X steps that preceded those. In other
	// words, we look for an apparent local convergence of the entropy measure of the MCMC.
//	std::cout << "\nstep   \tLogL       \tbest LogL\tMC step\n";
	newMeanL = -1e49;
	flagEq = false;
	bestL = d->getLikelihood();
	while (!flagEq) {
		oldMeanL = newMeanL;
		newMeanL = 0.0;
		for (int i = 0; i < 65536; i++) {
			if (!(d->monteCarloMove(dL, flagTaken))) {// Make a single MCMC move
				return false;
			}

			likeli = d->getLikelihood();			// get likelihood of this D
			if (likeli > bestL) {
				bestL = likeli;
			}	// store the current best likelihood
			newMeanL += likeli;

			if (t > 2147483640 || t < 0) {
				t = 1;
			} else {
				t++;
			}
		}
		d->refreshLikelihood();			// correct floating-point errors O(n)

		// Check if localized entropy appears to have stabilized; if you want to use a
		// different convergence criteria, this is where you would change things.
		if (std::abs(newMeanL - oldMeanL) / 65536.0 < 1.0
				&& (t > 10000 * ioparm.n || ioparm.flag_f_hrg)) {
			flagEq = true;
		}
	}
	return true;
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> bool UHRGPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::mcmcEquilibriumSample() {
	double dL, likeli, bestL;
	bool flagTaken;
	double ptest = 1.0 / 10.0; //(double)(4.0/(double)(ioparm.n));
	int thresh = 100 * ioparm.n;
	int t = 1;
	int sampleNum = 0;

	// Because moves in the dendrogram space are chosen (Monte Carlo) so that we sample dendrograms
	// with probability proportional to their likelihood, a likelihood-proportional sampling of
	// the dendrogram models would be equivalent to a uniform sampling of the walk itself. We would
	// still have to decide how often to sample the walk (at most once every n steps is recommended)
	// but for simplicity, the code here simply runs the MCMC itself. To actually compute something
	// over the set of sampled dendrogram models (in a Bayesian model averaging sense), you'll need
	// to code that yourself.

	bestL = d->getLikelihood();
	while (sampleNum < nbSamples) {
		for (int i = 0; i < 65536; i++) {
			// Make a single MCMC move
			if (!(d->monteCarloMove(dL, flagTaken))) {
				return false;
			}
			likeli = d->getLikelihood();				// get this likelihood
			if (likeli > bestL) {
				bestL = likeli;
			} // check if this logL beats best

			// We sample the dendrogram space every 1/ptest MCMC moves (on average).
			if (t > thresh && rng.getDouble(0, 1) < ptest) {
				sampleNum++;
				d->sampleAdjacencyLikelihoods();	// sample edge likelihoods
				if (sampleNum > nbSamples) {
					i = 65536;
				}
			}

			if (t > 2147483640 || t < 0) {
				t = 1;
			} else {
				t++;
			}
		}
		d->refreshLikelihood();			// correct floating-point errors O(n)
	}

	return true;
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> bool UHRGPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::toHRGNet() {

	// Add nodes
	ioparm.n = net->getNbNodes();
	g = new simpleGraph(ioparm.n, rng.getInt());// make new simpleGraph with n vertices
	d->g = new graph(ioparm.n);				// make new graph with n vertices
	d->g->setAdjacencyHistograms(nbBeans);		// setup adjacency histograms

	// Add edges to the graph
	ioparm.flag_timer = false;					// reset timer
	for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it) {
		auto a = Network::start(*it);
		auto b = Network::end(*it);
		if (!(g->addLink(a, b))) {
			throw std::runtime_error("Error adding edge");
		}
		if (!(g->addLink(b, a))) {
			throw std::runtime_error("Error adding edge");
		}
		if (!(d->g->addLink(a, b))) {
			throw std::runtime_error("Error adding edge");
		}
		if (!(d->g->addLink(b, a))) {
			throw std::runtime_error("Error adding edge");
		}
	}
	ioparm.m = g->getNumLinks();// store actual number of directional edges created
	ioparm.n = g->getNumNodes();			// store actual number of nodes used

	return true;
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UHRGPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::cleanup() {
	if (d != nullptr) {
		d->g->resetAllAdjacencies();			// clear A' for next trial
		d->g->resetLinks();					// clear G' for next trial
		d->resetDendrograph();				// clear D' for next trial
		delete d;
	}
	d = nullptr;
	delete g;
	g = nullptr;
	scores.clear();
}

#define UHRGPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef UHRGPREDICTOR_CPP


} /* namespace LinkPred */
