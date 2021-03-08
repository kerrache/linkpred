#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred::Simp;
int main() {
	// Create a prtedictor object
	Predictor p;
	// Load network from file
	p.loadnet("Zakarays_Karate_Club.edges"); 
	// Compute the score for the two edges (1, 34) and (26,34)
	std::vector<EdgeScore> esv = {{"1","34"},{"26","34"}};
	p.predADA(esv);
	// Print the scores
	for (auto it = esv.begin(); it != esv.end(); ++it) {
		std::cout << it->i << "\t" << it->j << "\t" << it->score << std::endl;
	}
	return 0;
}
