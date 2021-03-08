#include <linkpred.hpp>
#include <iostream>

using namespace LinkPred::Simp;

int main() {
	Predictor p; // Create a prtedictor object
	p.loadnet("Zakarays_Karate_Club.edges"); // Load network from file
	std::vector<EdgeScore> esv = { { "1", "34" }, { "26", "34" } };
	p.predKAB(esv);
	for (auto it = esv.begin(); it != esv.end(); ++it) {
		// Print the scores
		std::cout << it->i << "\t" << it->j << "\t" << it->score << std::endl;
	}
	return 0;
}
