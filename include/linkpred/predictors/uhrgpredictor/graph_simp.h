// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// graph_simp.h - graph data structure
// Copyright (C) 2006-2008 Aaron Clauset
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
// Simple graph data structure. The basic structure is an adjacency list of edges, along with degree
// information for the vertices.
// 
// ****************************************************************************************************

#ifndef GRAPH_SIMP_H
#define GRAPH_SIMP_H

#include <linkpred/predictors/uhrgpredictor/hrgds.h>
#include <linkpred/predictors/uhrgpredictor/rbtree.h>
#include "linkpred/utils/randomgen.hpp"
#include <stdio.h>
#include <string>
#include <stdlib.h>


// ******** Basic Structures ******************************************************************************

class simpleEdge {
public:
	int x;						// index of edge terminator
	simpleEdge* next;				// pointer to next elementd

	simpleEdge();
	~simpleEdge();
};

class simpleVert {
public:
	std::string name;					// (external) name of vertex
	int degree;					// degree of this vertex
	int group_true;				// index of vertex's true group

	simpleVert();
	~simpleVert();
};

class twoEdge {
public:
	int o;						// index of edge originator
	int x;						// index of edge terminator

	twoEdge();
	~twoEdge();
};

// ******** Graph Class with Edge Statistics *************************************************************

class simpleGraph {
public:
	simpleGraph(const int n, long int seed);
	~simpleGraph();

	bool addGroup(const int, const int);		// add group label to vertex i
	bool addLink(const int, const int);					// add (i,j) to graph
	bool doesLinkExist(const int, const int);// true if (i,j) is already in graph
	double getAdjacency(const int, const int);				// returns A(i,j)
	int getDegree(const int);					// returns degree of vertex i
	int getGroupLabel(const int);			// returns group label of vertex i
	std::string getName(const int);					// returns name of vertex i
	simpleEdge* getNeighborList(const int);		// returns edge list of vertex i
	simpleVert* getNode(const int);					// return pointer to a node
	int getNumGroups();								// returns num_groups
	int getNumLinks();									// returns m
	int getNumNodes();									// returns n
	bool setName(const int, const std::string);			// set name of vertex i

	void printAdjacencies();						// print adjacency matrix
	void printPairs();								// prints all edges in graph

private:
	simpleVert* nodes;			// list of nodes
	simpleEdge** nodeLink;			// linked list of neighbors to vertex
	simpleEdge** nodeLinkTail;		// pointers to tail of neighbor list
	double** A;				// adjacency matrix for this graph
	twoEdge* E;				// list of all edges (array)
	int n;				// number of vertices
	int m;				// number of directed edges
	int num_groups;		// number of bins in node histograms
	long int seed; // Seed
	LinkPred::RandomGen rng; // Mersenne Twister random number generator instance

	void QsortMain(block*, int, int);					// quicksort functions
	int QsortPartition(block*, int, int, int);
};

#endif
