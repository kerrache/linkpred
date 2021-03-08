#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred::Simp;
int main() {
	int nbRuns = 10;
	double edgeRemRatio = 0.1; 
	// Create an evaluator object
	Evaluator eval;
	// Add predictors to be evaluated
	eval.addCNE();
	eval.addADA();
	eval.addKAB();
	// Add performance measures
	eval.addROC();
	eval.addTPR();
	// Run experiment on the specified network
	eval.run("Zakarays_Karate_Club.edges", nbRuns, edgeRemRatio);
	// Print the header row
	auto res = eval.getPerfRes(0);
	for (auto it =  res.begin(); it != res.end(); ++it) {
		std::cout << it->name << "\t" ;
	}
	std::cout << "\n";
	// Print the results of each iteration
	for(int i = 0; i < nbRuns; i++) {
		auto res = eval.getPerfRes(i);
		for (auto it =  res.begin(); it != res.end(); ++it) {
			std::cout << it->res << "\t";
		}
		std::cout << "\n";
	}
	return 0;
}
