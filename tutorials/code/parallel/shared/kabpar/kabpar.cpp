#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main() {
	int k = 10; // Find top 10
	// Read network from file
	auto net = UNetwork<>::read("Infectious.edges"); 
	// Create an instance of the KAB predictor
	UKABPredictor<> p(net);
	// Enable parllelism
	p.setParallel(true);
	// Initialize predictor
	p.init();
	// Train predictor
	p.learn();
	// Allocate memory for storing scores
	std::vector<double> scores(k);
	// Create a vector to store edges
	std::vector<typename UNetwork<>::Edge> ev(k);
	// Predict top scores
	k = p.top(k, ev.begin(), scores.begin());
	for (int l = 0; l < k; l++) {
		auto i = net->start(ev[l]);
		auto j = net->end(ev[l]);
		std::cout << net->getLabel(i) << "\t" << net->getLabel(j) << "\t" << scores[l] << std::endl;
	}
	return 0;
}
