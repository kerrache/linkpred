#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred::Simp;
int main() {
	// We remove 10% of the edges
	double edgeRemRatio = 0.1;
	// We will not keep the network connected when removing edges
	bool keepConnected = false;
	// Seed of the random number generator
	long int seed = 0;
	// Create an Evaluator object
	Evaluator eval;
	// The ground truth network "Zakarays_Karate_Club.edges" is split into an observed network stored in "Zakarays_Train.edges" and a list of removed edges stored in "Zakarays_Test.edges"
	eval.genTestData("Zakarays_Karate_Club.edges", "Zakarays_Train.edges", "Zakarays_Test.edges", edgeRemRatio, keepConnected, seed);
	return 0;
}
