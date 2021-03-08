#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main() {
	int  k = 10;
	// Read network from file
	auto net = UNetwork<>::read("Zakarays_Karate_Club.edges"); 
	// Create a N2V encoder
	auto encoder =  std::make_shared<Node2Vec<>>(net, 777);
	// Create a logistic regresser
	auto classifier =  std::make_shared<LogisticRegresser<>>(0.001, 888);
	// Create an instance of the ECL predictor
	UECLPredictor<> p(net, encoder, classifier, 999);
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
