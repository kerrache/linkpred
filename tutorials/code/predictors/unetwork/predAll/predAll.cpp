#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main() {
	// Read network from file
	auto net = UNetwork<>::read("Zakarays_Karate_Club.edges"); 
	// Create an instance of the KAB predictor
	UKABPredictor<> p(net);
	// Initialize predictor
	p.init();
	// Train predictor
	p.learn();
	// Allocate memory for storing scores
	std::vector<double> scores(net->getNbNonEdges());
	// Predict the score of all non-existing edges
	auto its = p.predictNeg(scores.begin());
	// Print scores
	std::cout << "#Start\tEnd\tScore\n";
	int k = 0;
	for (auto it = its.first; it != its.second; ++it) {
		auto i = net->start(*it);
		auto j = net->end(*it);
		std::cout << net->getLabel(i) << "\t" << net->getLabel(j) << "\t" << scores[k++] << std::endl;
	}
	return 0;
}
