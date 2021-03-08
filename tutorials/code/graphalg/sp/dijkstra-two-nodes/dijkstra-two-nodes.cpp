#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main() {
	// Read network from file
	auto net = UNetwork<>::read("net-sp.edges");
	// Create an edge length (weight) map
	auto length = net->template createEdgeMapSP<double>();
	// Assign a length to every edge
	int i = 1;
	for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it, i++) {
		(*length)[*it] = (13 * i) % 3 + 1;
	}
	// Create a Dijkstra object
	Dijkstra<> dijkstra(net);
	// Register length map
	auto lengthMapId = dijkstra.registerLengthMap(length);
	// Find shortest path between node 1 and 6
	auto res = dijkstra.getShortestPath(net->getID("1"), net->getID("6"), lengthMapId);
	// Print path and distance
	auto path = res.first;
	auto dist = res.second;
	std::cout << "Path: ";
	for (auto it = path->begin(); it != path->end(); ++it) {
		std::cout << net->getLabel(*it) << " ";
	}
	std::cout << "\nDistance: " << dist << std::endl;
	return 0;
}
