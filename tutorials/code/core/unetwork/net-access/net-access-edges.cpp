#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main() {
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
	// Accessing edges
	std::cout << "Edges:\n";
	std::cout << "Start\tEnd\n";
	for (auto it = net.edgesBegin(); it != net.edgesEnd(); ++it) {
		std::cout << net.start(*it) << "\t" << net.end(*it) << std::endl;
	}
	// Neighbors
	for (std::size_t i = 0; i < net.getNbNodes(); i++) {
		std::cout << "Neighbors of " << i << std::endl;
		for (auto it = net.neighbBegin(i); it != net.neighbEnd(i); ++it) {
			std::cout << net.end(*it) << std::endl;
		}
	}
	return 0;
}
