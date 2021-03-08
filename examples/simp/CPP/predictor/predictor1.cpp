#include <linkpred.hpp>
#include <iostream>

using namespace LinkPred::Simp;

int main() {
	Predictor p; // Create a prtedictor object
	p.loadnet("Zakarays_Karate_Club.edges"); // Load network from file
	std::vector<EdgeScore> esv = p.predAllADA(); // Predict the score of all non-exisitng edges using Adamic Adar index.
	for (auto it = esv.begin(); it != esv.end(); ++it) {
		// Print the scores
		std::cout << it->i << "\t" << it->j << "\t" << it->score << std::endl;
	}
	return 0;
}
