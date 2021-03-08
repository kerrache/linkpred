#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main(int argc, char*argv[]) {
	// Initialize MPI
	MPI_Init(&argc, &argv);
	std::size_t k = 10;
	// Read network from file	
	auto net = UNetwork<>::read("Infectious.edges");
	// Create an RAL predictor
	URALPredictor<> p(net);
	// Enable distributed processing
	p.setDistributed(true);
	// Initialize predictor
	p.init();
	// Train predictor
	p.learn();
	// Allocate memory to store results
	std::vector<typename UNetwork<>::Edge> edges(k);
	std::vector<double> scores(k);
	// Find top k edges
	k = p.top(k, edges.begin(), scores.begin());
	int procID;
	// Get local process ID
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);
	// Print tyhe results
	if (procID == 0) {
		std::cout << "#Start\tEnd\tScore\n";
	}
	for (std::size_t i = 0; i < k; i++) {
		std::cout << net->getLabel(net->start(edges[i])) << "\t" << net->getLabel(net->end(edges[i])) << "\t" << scores[i] << std::endl;
	}
	// Finalize MPI
	MPI_Finalize();
	return 0;
}
