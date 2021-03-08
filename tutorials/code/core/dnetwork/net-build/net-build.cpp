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
	// Print the network
	net.print();
	return 0;
}
