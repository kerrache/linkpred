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
	// Create a node map that associates a double to every node
	auto nodeMap = net.template createNodeMap<double>();
	// Fill the map
	for (std::size_t i = 0; i < net.getNbNodes(); i++) {
		nodeMap[i] = i / 2.0;
	}
	// Access the map
	std::cout << "Label\tValue" << std::endl;
	for (std::size_t i = 0; i < net.getNbNodes(); i++) {
		std::cout << net.getLabel(i) << "\t" << nodeMap[i] << std::endl;
	}
	return 0;
}
