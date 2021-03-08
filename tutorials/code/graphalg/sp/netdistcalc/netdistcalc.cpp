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
	// Create a distance calculator. Here, we pass the option NetworkCache, that is we cache all distances
	ESPDistCalculator<> calc(dijkstra, length, NetworkCache);
	// Print all distances
	std::cout << "i\tj\tDist" << std::endl;
	for (unsigned int i = 0; i < net->getNbNodes(); i++) {
		for (unsigned int j = 0; j < i; j++) {
			std::cout << net->getLabel(i) << "\t" << net->getLabel(j) << "\t" << calc.getDist(i, j).first << std::endl;
		}
	}
	return 0;
}
