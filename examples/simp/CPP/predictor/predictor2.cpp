#include <linkpred.hpp>
#include <iostream>

using namespace LinkPred::Simp;

int main() {
	Predictor p;
	p.loadnet("Zakarays_Karate_Club.edges");
	auto esv = p.predAllKAB();
	for (auto it = esv.begin(); it != esv.end(); ++it) {
		std::cout << it->i << "\t" << it->j << "\t" << it->score << std::endl;
	}

	return 0;
}
