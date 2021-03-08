#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main(int argc, char*argv[]) {
	// Create a network object 
	UNetwork<> net;
	// Add nodes
	net.addNode("A");
	net.addNode("B");
	net.addNode("C");
	net.addNode("D");
	// Add edges
	net.addEdge(net.getID("A"), net.getID("B"));
	net.addEdge(net.getID("B"), net.getID("C"));
	net.addEdge(net.getID("C"), net.getID("D"));
	net.addEdge(net.getID("D"), net.getID("A"));
	// Assemble the network
	net.assemble();
	// Create a node map that associates an integer to every edge
	auto edgeMap = net.template createEdgeMap<int>();
	edgeMap[net.makeEdge(0, 1)] = 3;
	edgeMap[net.makeEdge(1, 2)] = 2;
	edgeMap[net.makeEdge(2, 3)] = 5;
	edgeMap[net.makeEdge(3, 0)] = 4;
	// Access the map
	std::cout << "Start\tEnd\tValue" << std::endl;
	for (auto it = net.edgesBegin(); it != net.edgesEnd(); ++it) {
		auto i = net.start(*it);
		auto j = net.end(*it);
		std::cout << net.getLabel(i) << "\t" << net.getLabel(j) << "\t" << edgeMap.at(*it) << std::endl;
	}
	return 0;
}
