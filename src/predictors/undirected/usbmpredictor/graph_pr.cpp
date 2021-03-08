/*
 * graph_pr.cpp
 *
 *  Created on: Sep 6, 2016
 *      Author: Said Kerrache
 */

#include "linkpred/predictors/undirected/uhrgpredictor/graph_pr.h"
#include <cmath>

edge::edge() {
	x = -1;
	next = nullptr;
	h = nullptr;
	total_weight = 0.0;
	obs_count = 0;
}
edge::~edge() {
	if (h != nullptr) {
		delete[] h;
	}
	h = nullptr;
}
void edge::setHistogram(const int k) {
	h = new double[k];
	for (int i = 0; i < k; i++) {
		h[i] = 0.0;
	}
}
void edge::resetHistogram(const int k) {
	total_weight = 0.0;
	obs_count = 0;
	for (int i = 0; i < k; i++) {
		h[i] = 0.0;
	}
}

vert::vert() {
	name = "";
	degree = 0;
}
vert::~vert() {
}

graph::graph(const int size) {
	n = size;
	m = 0;
	num_bins = 0;
	bin_resolution = 0.0;
	nodes = new vert[n];
	nodeLink = new edge*[n];
	nodeLinkTail = new edge*[n];
	A = new double**[n];
	for (int i = 0; i < n; i++) {
		nodeLink[i] = nullptr;
		nodeLinkTail[i] = nullptr;
		A[i] = new double*[n];
	}
	obs_count = 0;
	total_weight = 0.0;
}

graph::~graph() {
	edge *curr, *prev;
	for (int i = 0; i < n; i++) {
		curr = nodeLink[i];
		while (curr != nullptr) {
			prev = curr;
			curr = curr->next;
			delete prev;						// deletes edge histogram, too
		}
		for (int j = 0; j < n; j++) {
			delete[] A[i][j];
		}
		delete[] A[i];
	}
	delete[] A;
	A = nullptr;
	delete[] nodeLink;
	nodeLink = nullptr;
	delete[] nodeLinkTail;
	nodeLinkTail = nullptr;
	delete[] nodes;
	nodes = nullptr;	// deletes node histogram, too
}

// ********************************************************************************************************

bool graph::addLink(const int i, const int j) {
	// Adds the directed edge (i,j) to the adjacency list for v_i
	edge* newedge;
	if (i >= 0 && i < n && j >= 0 && j < n) {
		newedge = new edge;
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
		return true;
	} else {
		return false;
	}
}

// ********************************************************************************************************

bool graph::addAdjacencyObs(const int i, const int j, const double probability,
		const double size) {
	// Adds the observation obs to the histogram of the edge (i,j)
	// Note: user must manually add observation to edge (j,i) by calling this function with that argument
//	char pauseme;
	if (bin_resolution > 0.0 && probability >= 0.0 && probability <= 1.0
			&& size >= 0.0 && size <= 1.0 && i >= 0 && i < n && j >= 0
			&& j < n) {
		int index = (int) (std::round(probability / bin_resolution));
		if (index < 0) {
			index = 0;
		} else if (index > num_bins) {
			index = num_bins;
		}
		// Add the weight to the proper probability bin
		if (A[i][j][index] < 0.5) {
			A[i][j][index] = 1.0;
		} else {
			A[i][j][index] += 1.0;
		}
		return true;
	}
	return false;
}

// ********************************************************************************************************

void graph::addAdjacencyEnd() {
	// We need to also keep a running total of how much weight has been added
	// to the histogram, and the number of observations in the histogram.
	if (obs_count == 0) {
		total_weight = 1.0;
		obs_count = 1;
	} else {
		total_weight += 1.0;
		obs_count++;
	}
	return;
}

// ********************************************************************************************************

bool graph::doesLinkExist(const int i, const int j) {
	// This function determines if the edge (i,j) already exists in the adjacency list of v_i
	edge* curr;
	if (i >= 0 && i < n && j >= 0 && j < n) {
		curr = nodeLink[i];
		while (curr != nullptr) {
			if (curr->x == j) {
				return true;
			}
			curr = curr->next;
		}
	}
	return false;
}

// ********************************************************************************************************

int graph::getDegree(const int i) {
	if (i >= 0 && i < n) {
		return nodes[i].degree;
	} else {
		return -1;
	}
}
std::string graph::getName(const int i) {
	if (i >= 0 && i < n) {
		return nodes[i].name;
	} else {
		return "";
	}
}
// NOTE: The following three functions return addresses; deallocation of returned object is dangerous
edge* graph::getNeighborList(const int i) {
	if (i >= 0 && i < n) {
		return nodeLink[i];
	} else {
		return nullptr;
	}
}
double* graph::getAdjacencyHist(const int i, const int j) {
	if (i >= 0 && i < n && j >= 0 && j < n) {
		return A[i][j];
	} else {
		return nullptr;
	}
}
// END-NOTE

// ********************************************************************************************************

double graph::getAdjacencyAverage(const int i, const int j) {
	double average = 0.0;
//	double temp[num_bins];
//	double tot = 0.0;
	if (i != j) {
		for (int k = 0; k < num_bins; k++) {
			if (A[i][j][k] > 0.0) {
				average += (A[i][j][k] / total_weight)
						* ((double) (k) * bin_resolution);
			}
		}
	}
	return average;
}

// ********************************************************************************************************

double graph::getBinResolution() {
	return bin_resolution;
}
int graph::getNbBeans() {
	return num_bins;
}
int graph::getNumLinks() {
	return m;
}
int graph::getNumNodes() {
	return n;
}
double graph::getTotalWeight() {
	return total_weight;
}

// ********************************************************************************************************

void graph::printPairs() {
	edge* curr;
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
	std::cout << edgeCount << " edges total.\n";
	return;
}

// ********************************************************************************************************

void graph::printAdjacencies() {
	double average;
	for (int i = 0; i < n; i++) {
		std::cout << "[" << i << "]";
		for (int j = 0; j < n; j++) {
			average = 0.0;
			for (int k = 0; k < num_bins; k++) {
				if (A[i][j][k] > 0.0) {
					average += (A[i][j][k] / total_weight)
							* ((double) (k) * bin_resolution);
				}
			}
			average = round(average * 100.0) / 100.0;
			std::cout << " " << average;
		}
		std::cout << std::endl;
	}
	return;
}

// ********************************************************************************************************

void graph::printAdjacencyHists() {
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			std::cout << "(" << i << " " << j << ")\t[ ";
			for (int k = 0; k < num_bins; k++) {
				std::cout << A[i][j][k] << " ";
			}
			std::cout << "]\n";
		}
	}
	return;
}

void graph::printAdjacencyHist(const int i, const int j) {
	double average = 0.0;
	std::cout << "(" << i << " " << j << ")\t[ ";
	for (int k = 0; k < num_bins; k++) {
		if (A[i][j][k] > 0.0) {
			average += (A[i][j][k] / total_weight)
					* ((double) (k) * bin_resolution);
		}
		std::cout << A[i][j][k] << " ";
	}
	std::cout << "] tot = ";
	std::cout << total_weight << "\tmean = " << average << "\n";
	return;
}

// ********************************************************************************************************

void graph::resetAllAdjacencies() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < num_bins; k++) {
				A[i][j][k] = 0.0;
			}
		}
	}
	obs_count = 0;
	total_weight = 0.0;
	return;
}

// ********************************************************************************************************

void graph::resetAdjacencyHistogram(const int i, const int j) {
	if (i >= 0 && i < n && j >= 0 && j < n) {
		for (int k = 0; k < num_bins; k++) {
			A[i][j][k] = 0.0;
		}
	}
	return;
}

// ********************************************************************************************************

void graph::resetLinks() {
	edge *curr, *prev;
	for (int i = 0; i < n; i++) {
		curr = nodeLink[i];
		while (curr != nullptr) {
			prev = curr;
			curr = curr->next;
			delete prev;
		}
		nodeLink[i] = nullptr;
		nodeLinkTail[i] = nullptr;
		nodes[i].degree = 0;
	}
	m = 0;
	return;
}

// ********************************************************************************************************

void graph::setAdjacencyHistograms(const int bin_count) {
	// For all possible adjacencies, setup an edge histograms
	num_bins = bin_count + 1;
	bin_resolution = 1.0 / (double) (bin_count);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = new double[num_bins];
			for (int k = 0; k < num_bins; k++) {
				A[i][j][k] = 0.0;
			}
		}
	}
	return;
}

// ********************************************************************************************************

bool graph::setName(const int i, const std::string text) {
	if (i >= 0 && i < n) {
		nodes[i].name = text;
		return true;
	} else {
		return false;
	}
}

