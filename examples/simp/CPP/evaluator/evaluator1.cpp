#include <linkpred.hpp>
#include <iostream>

using namespace LinkPred::Simp;

int main() {
	int nbRuns = 10;
	Evaluator eval;
	eval.addCNE();
	eval.addADA();
	eval.addKAB("KAB-2", 2); // Horizon limit = 2
	eval.addKAB("KAB-4", 4); // Horizon limit = 4
	eval.addROC();
	eval.addTPR();
	eval.run("Zakarays_Karate_Club.edges", nbRuns);
	// Print the header row
	auto res = eval.getPerfRes(0);
	for (auto it = res.begin(); it != res.end(); ++it) {
		std::cout << it->name << "\t";
	}
	std::cout << "\n";
	// Print the results of each iteration
	for (int i = 0; i < nbRuns; i++) {
		auto res = eval.getPerfRes(i);
		for (auto it = res.begin(); it != res.end(); ++it) {
			std::cout << it->res << "\t";
		}
		std::cout << "\n";
	}
	return 0;
}
