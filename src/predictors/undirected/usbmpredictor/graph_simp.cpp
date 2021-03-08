/*
 * graph_simp.cpp
 *
 *  Created on: Sep 6, 2016
 *      Author: Said Kerrache
 */

#include "linkpred/predictors/undirected/uhrgpredictor/graph_simp.h"

simpleEdge::simpleEdge() {
	x = -1;
	next = nullptr;
}
simpleEdge::~simpleEdge() {
}

simpleVert::simpleVert() {
	name = "";
	degree = 0;
	group_true = -1;
}
simpleVert::~simpleVert() {
}

twoEdge::twoEdge() {
	o = -1;
	x = -1;
}
twoEdge::~twoEdge() {
}

simpleGraph::simpleGraph(const int size, long int seed) :
		seed(seed), rng(seed) {
	n = size;
	m = 0;
	num_groups = 0;
	nodes = new simpleVert[n];
	nodeLink = new simpleEdge*[n];
	nodeLinkTail = new simpleEdge*[n];
	A = new double*[n];
	for (int i = 0; i < n; i++) {
		nodeLink[i] = nullptr;
		nodeLinkTail[i] = nullptr;
		A[i] = new double[n];
		for (int j = 0; j < n; j++) {
			A[i][j] = 0.0;
		}
	}
	E = nullptr;
}

simpleGraph::~simpleGraph() {
	simpleEdge *curr, *prev;
	for (int i = 0; i < n; i++) {
		curr = nodeLink[i];
		delete[] A[i];
		while (curr != nullptr) {
			prev = curr;
			curr = curr->next;
			delete prev;
		}
	}
	curr = nullptr;
	prev = nullptr;
	if (E != nullptr) {
		delete[] E;
		E = nullptr;
	}
	delete[] A;
	A = nullptr;
	delete[] nodeLink;
	nodeLink = nullptr;
	delete[] nodeLinkTail;
	nodeLinkTail = nullptr;
	delete[] nodes;
	nodes = nullptr;
}

// ********************************************************************************************************

bool simpleGraph::addGroup(const int i, const int group_index) {
	if (i >= 0 && i < n) {
		nodes[i].group_true = group_index;
		return true;
	} else {
		return false;
	}
}

// ********************************************************************************************************

bool simpleGraph::addLink(const int i, const int j) {
	// Adds the directed edge (i,j) to the adjacency list for v_i
	simpleEdge* newedge;
	if (i >= 0 && i < n && j >= 0 && j < n) {
		A[i][j] = 1.0;
		newedge = new simpleEdge;
		newedge->x = j;
		if (nodeLink[i] == nullptr) {			// first neighbor
			nodeLink[i] = newedge;
			nodeLinkTail[i] = newedge;
			nodes[i].degree = 1;
		} else {							// subsequent neighbor
			nodeLinkTail[i]->next = newedge;
			nodeLinkTail[i] = newedge;
			nodes[i].degree++;
		}
		m++;								// increment edge count
		newedge = nullptr;
		return true;
	} else {
		return false;
	}
}

// ********************************************************************************************************

bool simpleGraph::doesLinkExist(const int i, const int j) {
	// This function determines if the edge (i,j) already exists in the adjacency list of v_i
	if (i >= 0 && i < n && j >= 0 && j < n) {
		if (A[i][j] > 0.1) {
			return true;
		} else {
			return false;
		}
	} else {
		return false;
	}
	return false;
}

// ********************************************************************************************************

double simpleGraph::getAdjacency(const int i, const int j) {
	if (i >= 0 && i < n && j >= 0 && j < n) {
		return A[i][j];
	} else {
		return -1.0;
	}
}
int simpleGraph::getDegree(const int i) {
	if (i >= 0 && i < n) {
		return nodes[i].degree;
	} else {
		return -1;
	}
}
int simpleGraph::getGroupLabel(const int i) {
	if (i >= 0 && i < n) {
		return nodes[i].group_true;
	} else {
		return -1;
	}
}
std::string simpleGraph::getName(const int i) {
	if (i >= 0 && i < n) {
		return nodes[i].name;
	} else {
		return "";
	}
}
// NOTE: The following three functions return addresses; deallocation of returned object is dangerous
simpleEdge* simpleGraph::getNeighborList(const int i) {
	if (i >= 0 && i < n) {
		return nodeLink[i];
	} else {
		return nullptr;
	}
}
// END-NOTE

// ********************************************************************************************************

int simpleGraph::getNumGroups() {
	return num_groups;
}
int simpleGraph::getNumLinks() {
	return m;
}
int simpleGraph::getNumNodes() {
	return n;
}
simpleVert* simpleGraph::getNode(const int i) {
	if (i >= 0 && i < n) {
		return &nodes[i];
	} else {
		return nullptr;
	}
}

// ********************************************************************************************************

void simpleGraph::printAdjacencies() {
	for (int i = 0; i < n; i++) {
		std::cout << "[" << i << "]";
		for (int j = 0; j < n; j++) {
			std::cout << " " << A[i][j];
		}
		std::cout << std::endl;
	}
	return;
}

// ********************************************************************************************************

void simpleGraph::printPairs() {
	simpleEdge* curr;
	int edgeCount = 0;
	for (int i = 0; i < n; i++) {
		std::cout << "[" << i << "]\t";
		curr = nodeLink[i];
		while (curr != nullptr) {
			std::cout << curr->x << "\t";
			edgeCount++;
			curr = curr->next;
		}
		std::cout << "\n";
	}
	curr = nullptr;
	std::cout << edgeCount << " edges total.\n";
	return;
}

// ********************************************************************************************************

bool simpleGraph::setName(const int i, const std::string text) {
	if (i >= 0 && i < n) {
		nodes[i].name = text;
		return true;
	} else {
		return false;
	}
}

// ********************************************************************************************************

void simpleGraph::QsortMain(block* array, int left, int right) {
	if (right > left) {
		int pivot = left;
		int part = QsortPartition(array, left, right, pivot);
		QsortMain(array, left, part - 1);
		QsortMain(array, part + 1, right);
	}
	return;
}

int simpleGraph::QsortPartition(block* array, int left, int right, int index) {
	block p_value, temp;
	p_value.x = array[index].x;
	p_value.y = array[index].y;

	// swap(array[p_value], array[right])
	temp.x = array[right].x;
	temp.y = array[right].y;
	array[right].x = array[index].x;
	array[right].y = array[index].y;
	array[index].x = temp.x;
	array[index].y = temp.y;

	int stored = left;
	for (int i = left; i < right; i++) {
		if (array[i].x <= p_value.x) {
			// swap(array[stored], array[i])
			temp.x = array[i].x;
			temp.y = array[i].y;
			array[i].x = array[stored].x;
			array[i].y = array[stored].y;
			array[stored].x = temp.x;
			array[stored].y = temp.y;
			stored++;
		}
	}
	// swap(array[right], array[stored])
	temp.x = array[stored].x;
	temp.y = array[stored].y;
	array[stored].x = array[right].x;
	array[stored].y = array[right].y;
	array[right].x = temp.x;
	array[right].y = temp.y;

	return stored;
}

