#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main() {
	// Read network from file
	auto net = UNetwork<>::read("net-traversal.edges");
	// Create a BFS objec
	BFS<> bfs(net);
	// We collect nodes during traversal
	Collector<> col;
	// We start traversal at node 1
	bfs.traverse(net->getID("1"), col);
	// Retrieve the set of visited nodes
	auto visited = col.getVisited();
	// Print visited nodes
	std::cout << "Visited nodes:" << std::endl;
	while (!visited.empty()) {
		auto i = visited.front();
		visited.pop();
		std::cout << net->getLabel(i) << std::endl;
	}
	return 0;
}
