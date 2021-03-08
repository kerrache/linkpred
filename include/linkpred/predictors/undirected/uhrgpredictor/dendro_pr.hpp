// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// dendro_pr.h - hierarchical random graph (hrg) data structure
// Copyright (C) 2005-2008 Aaron Clauset
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
// See http://www.gnu.org/licenses/gpl.txt for more details.
// 
// ****************************************************************************************************
// Author       : Aaron Clauset  ( aaronc@santafe.edu | http://www.santafe.edu/~aaronc/ )
// Collaborators: Cristopher Moore and Mark E.J. Newman
// Project      : Hierarchical Random Graphs
// Location     : University of New Mexico, Dept. of Computer Science AND Santa Fe Institute
// Created      : 21 June 2006
// Modified     : 23 December 2007 (cleaned up for public consumption)
//
// ****************************************************************************************************
// 
// Maximum likelihood dendrogram data structure. This is the heart of the HRG algorithm: all
// manipulations are done here and all data is stored here. The data structure uses the separate
// graph data structure to store the basic adjacency information (in a dangerously mutable way).
// 
// ****************************************************************************************************

#ifndef DENDRO_PR_H
#define DENDRO_PR_H

#include <linkpred/predictors/undirected/uhrgpredictor/graph_pr.h>
#include <linkpred/predictors/undirected/uhrgpredictor/hrgds.h>
#include <linkpred/predictors/undirected/uhrgpredictor/rbtree.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>

#include "linkpred/utils/randomgen.hpp"

// ********************************************************************************************************
// ******** Basic Structures ******************************************************************************

enum {
	DENDRO, GRAPH, LEFT, RIGHT
};

// ********************************************************************************************************
// ******** Internal Edge Class ***************************************************************************
// The usefulness of this data structure is to provide an easy to way maintain the set of internal edges,
// in the dendrogram D. It allows for the selection of a random internal edge in O(1) time, and it takes 
// O(1) time to update its structure given an internal move.

class interns {
private:
	ipair* edgelist;					// list of internal edges represented
	int** indexLUT;			// table of indices of internal edges in edgelist
	int q;						// number of internal edges
	int count;			// (for adding edges) edgelist index of new edge to add
	long int seed; // Seed
	LinkPred::RandomGen rng; // Mersenne Twister random number generator instance

public:
	interns(const int n, long int seed);
	~interns();

	bool addEdge(const int, const int, const short int); // add an internal edge, O(1)
	ipair* getEdge(const int);		// returns the ith edge of edgelist, O(1)
	ipair* getRandomEdge();	// returns a uniformly random internal edge, O(1)
	void printEdgeList();						// writes edgelist to terminal
	bool swapEdges(const int, const int, const short int, const int, const int,
			const short int);
	// swaps two edges, O(1)
};

// ********************************************************************************************************

// ********************************************************************************************************
// ******** Tree elementd Class ***************************************************************************

class elementd {
public:
	short int type;				// either DENDRO or GRAPH
	double logL;			// log-likelihood contribution of this internal node
	double p;	// probability p_i that an edge exists between L and R subtrees
	int e;					// number of edges between L and R subtrees
	int n;					// number of leafs in subtree rooted here
	int label;				// subtree label: smallest leaf index
	int index;				// index in containing array

	elementd *M;					// pointer to parent node
	elementd *L;					// pointer for L subtree
	elementd *R;					// pointer for R subtree

	elementd();
	~elementd();
};

// ********************************************************************************************************
// ******** Dendrogram Class ******************************************************************************

class dendro {

private:
	elementd* root = nullptr;			// root of the dendrogram
	elementd* internal = nullptr;// array of n-1 internal vertices (the dendrogram D)
	elementd* leaf = nullptr;		// array of n   leaf vertices (the graph G)
	int n = 0;				// number of leaf vertices to allocate
	interns* d = nullptr;			// list of internal edges of dendrogram D
	list** paths = nullptr;			// array of path-lists from root to leaf
	double L = 0;				// log-likelihood of graph G given dendrogram D
	long int seed; // Seed
	LinkPred::RandomGen rng; // Mersenne Twister random number generator instance
	rbtree subtreeL, subtreeR;	// trees for computeEdgeCount() function

	void binarySearchInsert(elementd*, elementd*);// insert node i according to binary search property
	list* binarySearchFind(const double);		// return path to root from leaf
	int computeEdgeCount(const int, const short int, const int,
			const short int); // compute number of edges between two internal subtrees
	elementd* findCommonAncestor(list**, const int, const int);	// find internal node of D that is common ancestor of i,j
	void printSubTree(elementd*);			// display the subtree rooted at z
	list* reversePathToRoot(const int);	// return reverse of path to leaf from root
	void QsortMain(block*, int, int);					// quicksort functions
	int QsortPartition(block*, int, int, int);

public:
	graph* g = nullptr;					// underlying G (dangerously accessible)

	dendro(long int seed) :
			seed(seed), rng(seed) {
	}
	~dendro();								// constructor / destructor
	void buildDendrogram();						// build dendrogram from g
	void clearDendrograph();// delete dendrograph in prep for importDendrogramStructure
	double getLikelihood();					// return likelihood of G given D
	bool importDendrogramStructure(const std::string);// read dendrogram structure from file
	bool monteCarloMove(double&, bool&);			// make single MCMC move
	void refreshLikelihood();			// force refresh of log-likelihood value
	void sampleAdjacencyLikelihoods();// sample dendrogram edge likelihoods and update edge histograms
	void resetDendrograph();				// reset the dendrograph structures
	void printDendrogram();			// write dendrogram structure to terminal
};

#endif
