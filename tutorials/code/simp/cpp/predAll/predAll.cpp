#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred::Simp;
int main() {
	// Create a prtedictor object
	Predictor p;
	// Load network from file
	p.loadnet("Zakarays_Karate_Club.edges"); 
	// Predict the score of all non-existing edges using Adamic Adar index
	std::vector<EdgeScore> esv = p.predAllADA(); 
	// Print the scores
	for (auto it = esv.begin(); it != esv.end(); ++it) {
		std::cout << it->i << "\t" << it->j << "\t" << it->score << std::endl;
	}
	return 0;
}
