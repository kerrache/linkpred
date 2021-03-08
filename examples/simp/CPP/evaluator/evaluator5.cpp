#include <linkpred.hpp>
#include <iostream>

using namespace LinkPred::Simp;

int main() {
	double edgeRemRatio = 0.1;
	bool keepConnected = false;
	long int seed = 0;
	Evaluator eval;
	eval.genTestData("Zakarays_Karate_Club.edges", "Zakarays_Train.edges",
			"Zakarays_Test.edges", edgeRemRatio, keepConnected, seed);
	return 0;
}
