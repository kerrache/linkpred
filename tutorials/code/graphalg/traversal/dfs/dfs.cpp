#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main() {
	// Read network from file
	auto net = UNetwork<>::read("net-traversal.edges");
	// Create a DFS objec
	DFS<UNetwork<>, Counter<>> dfs(net);
	// We count nodes during traversal
	Counter<> counter;
	// We start traversal at node 1
	dfs.traverse(net->getID("1"), counter);
	// Print the number of visited nodes
	std::cout << "DFS visited " << counter.getCount() << " nodes" << std::endl;
	return 0;
}
