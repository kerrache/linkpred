#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred::Simp;
int main() {
	int k = 10;
	// Create a prtedictor object
	Predictor p;
	// Load network from file
	p.loadnet("Zakarays_Karate_Club.edges");
	// Predict top k scores using an encoder-classifier predictor with Node2Vec as encoder and logistic regression as classifier
	auto esv = p.predTopECL(k, "N2V", "LGR");
	// Print scores
	std::cout << "N2V-LGR\n";
	for (auto it = esv.begin(); it != esv.end(); ++it) {
		std::cout << it->i << "\t" << it->j << "\t" << it->score << std::endl;
	}
	// Predict top k scores using an encoder-classifier predictor with LINE as encoder and feed-forward neural network as classifier
	esv = p.predTopECL(k, "LIN", "FFN"); // FFN requires mlpack
	// Print scores
	std::cout << "LIN-FFN\n";
	for (auto it = esv.begin(); it != esv.end(); ++it) {
		std::cout << it->i << "\t" << it->j << "\t" << it->score << std::endl;
	}
	return 0;
}
