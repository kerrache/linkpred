#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred::Simp;
int main() {
	int k = 10;
	// Create a prtedictor object
	Predictor p;
	// Load network from file
	p.loadnet("Zakarays_Karate_Club.edges");
	// Predict top k scores using an encoder-classifier predictor with Node2Vec as encoder and L2 similarity
	auto esv = p.predTopESM(k, "N2V", "L2");
	// Print scores
	std::cout << "N2V-L2\n";
	for (auto it = esv.begin(); it != esv.end(); ++it) {
		std::cout << it->i << "\t" << it->j << "\t" << it->score << std::endl;
	}
	// Predict top k scores using an encoder-similarity measure predictor with LINE as encoder and cosine similarity
	esv = p.predTopESM(k, "LIN", "CSM");
	// Print scores
	std::cout << "LIN-CSM\n";
	for (auto it = esv.begin(); it != esv.end(); ++it) {
		std::cout << it->i << "\t" << it->j << "\t" << it->score << std::endl;
	}
	return 0;
}
