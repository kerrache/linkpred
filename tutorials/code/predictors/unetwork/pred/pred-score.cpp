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
	// Compute the score for the two edges (1, 34) and (26,34)
	double sc = p.score(net->makeEdge(net->getID("1"), net->getID("34"))); 
	std::cout << "1\t34\t" << sc << std::endl;
	sc = p.score(net->makeEdge(net->getID("26"), net->getID("34"))); 
	std::cout << "26\t34\t" << sc << std::endl;	
	return 0;
}
