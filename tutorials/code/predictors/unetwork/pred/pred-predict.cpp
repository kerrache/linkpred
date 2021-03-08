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
	// Create a vector to store edges
	std::vector<typename UNetwork<>::Edge> ev;
	// Push the two edges (1, 34) and (26, 34)
	ev.push_back(net->makeEdge(net->getID("1"), net->getID("34"))); 
	ev.push_back(net->makeEdge(net->getID("26"), net->getID("34"))); 
	// Allocate memory for storing scores
	std::vector<double> scores(2);
	// Predict the scores
	p.predict(ev.begin(), ev.end(), scores.begin());
	// Print scores
	std::cout << "#Start\tEnd\tScore\n";
	int k = 0;
	for (auto it = ev.begin(); it != ev.end(); ++it) {
		auto i = net->start(*it);
		auto j = net->end(*it);
		std::cout << net->getLabel(i) << "\t" << net->getLabel(j) << "\t" << scores[k++] << std::endl;
	}

	return 0;
}
