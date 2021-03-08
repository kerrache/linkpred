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
	// Compute the distance from node 1 to all other nodes
	auto distMap = dijkstra.getDist(net->getID("1"), lengthMapId);
	// Print distances
	std::cout << "Target\tDist\tNumber of nodes in the path" << std::endl;
	for (auto it = net->nodesBegin(); it != net->nodesEnd(); ++it) {
		auto res = distMap->at(it->first);
		std::cout << it->second << "\t" << res.first << "\t" << res.second << std::endl;
	}
	return 0;
}
