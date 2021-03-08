#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main() {
	// Create a network object 
	DNetwork<> net;
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
	// Accessing nodes
	std::cout << "Nodes:\n";
	std::cout << "ID\tLabel\n";
	for (auto it = net.nodesBegin(); it != net.nodesEnd(); ++it) {
		std::cout << it->first << "\t" << it->second << std::endl;
	}
	std::cout << "Translating labels to IDs\n";
	std::cout << "A -> " << net.getID("A") << std::endl;
	std::cout << "B -> " << net.getID("B") << std::endl;
	std::cout << "C -> " << net.getID("C") << std::endl;
	std::cout << "D -> " << net.getID("D") << std::endl;
	std::cout << "Translating IDs to labels\n";
	for (std::size_t i = 0; i < net.getNbNodes(); i++) {
		std::cout << i << " -> " << net.getLabel(i) << std::endl;
	}	
	return 0;
}
