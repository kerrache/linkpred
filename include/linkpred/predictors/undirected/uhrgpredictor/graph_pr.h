// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// graph_pr.h - graph data structure
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
// Graph data structure for maximum likelihood hrgs. The basic structure is an adjacency list of
// edges; however, many additional pieces of metadata are stored as well. Each node stores its
// external name and its degree. Each edge stores a histogram of the probabilities assigned to it 
// by the dendrogram structure. Generally, edges are directional, an adjacency (i,j) should also be 
// stored as (j,i) in order to maintain the undirectedness of the graph.
// 
// ****************************************************************************************************

#ifndef GRAPH_PR_H
#define GRAPH_PR_H

#include <linkpred/predictors/undirected/uhrgpredictor/rbtree.h>
#include <stdio.h>
#include <string>
#include <stdlib.h>

// ******** Basic Structures ******************************************************************************

class edge {
public:
	int x;						// index of edge terminator
	double* h;						// (histogram) weights of edge existence
	double total_weight;				// (histogram) total weight observed
	int obs_count;			// (histogram) number of observations in histogram
	edge* next;					// pointer to next elementd

	edge();
	~edge();
	void setHistogram(const int);		// allocate / setup histogram
	void resetHistogram(const int);	// clear histogram data
};

class vert {
public:
	std::string name;					// (external) name of vertex
	int degree;					// degree of this vertex

	vert();
	~vert();
};

// ******** Graph Class with Edge Statistics *************************************************************

class graph {
public:
	graph(const int);
	~graph();

	bool addLink(const int, const int);					// add (i,j) to graph
	bool addAdjacencyObs(const int, const int, const double, const double);	// add weight to (i,j)'s histogram
	void addAdjacencyEnd();					// add to obs_count and total_weight
	bool doesLinkExist(const int, const int);// true if (i,j) is already in graph
	int getDegree(const int);					// returns degree of vertex i
	std::string getName(const int);					// returns name of vertex i
	edge* getNeighborList(const int);			// returns edge list of vertex i
	double* getAdjacencyHist(const int, const int);	// return ptr to histogram of edge (i,j)
	double getAdjacencyAverage(const int, const int);// return average value of adjacency A(i,j)
	double getBinResolution();						// returns bin_resolution
	int getNbBeans();									// returns num_bins
	int getNumLinks();									// returns m
	int getNumNodes();									// returns n
	double getTotalWeight();							// returns total_weight
	void resetAdjacencyHistogram(const int, const int);	// reset edge (i,j)'s histogram
	void resetAllAdjacencies();						// reset all edge histograms
	void resetLinks();							// clear all links from graph
	void setAdjacencyHistograms(const int);			// allocate edge histograms
	bool setName(const int, const std::string);			// set name of vertex i

	void printPairs();								// prints all edges in graph
	void printAdjacencies();				// print average adjacency values
	void printAdjacencyHist(const int i, const int j);// print one edge histogram
	void printAdjacencyHists();			// prints edges with their histograms

private:
	vert* nodes;			// list of nodes
	edge** nodeLink;			// linked list of neighbors to vertex
	edge** nodeLinkTail;		// pointers to tail of neighbor list
	double*** A;				// stochastic adjacency matrix for this graph
	int obs_count;		// number of observations in A
	double total_weight;		// total weight added to A
	int n;				// number of vertices
	int m;				// number of directed edges
	int num_bins;			// number of bins in edge histograms
	double bin_resolution;	// width of histogram bin
};

#endif
