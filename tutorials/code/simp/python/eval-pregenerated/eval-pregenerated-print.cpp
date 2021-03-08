#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred::Simp;
int main() {
	// Create an evaluator object
	Evaluator eval;
	// Add predictors to be evaluated
	eval.addADA();
	eval.addRAL();
	// Add performance measures
	eval.addPR();
	eval.addTPR();
	// Run experiment on the specified network
	eval.run("Zakarays_Karate_Club_Train.edges", "Zakarays_Karate_Club_Test.edges");
	auto res = eval.getPerfRes(0);
	for (auto it =  res.begin(); it != res.end(); ++it) {
		std::cout << it->name << "\t" ;
	}
	std::cout << "\n";
	for (auto it =  res.begin(); it != res.end(); ++it) {
		std::cout << it->res << "\t" ;
	}
	std::cout << "\n";
	return 0;
}
