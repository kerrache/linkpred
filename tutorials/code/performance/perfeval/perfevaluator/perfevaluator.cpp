#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main() {
	// Remove 10% of the edges
	double remRatio = 0.1;
	long int seed = 888;
	// Read network from file
	auto refNet = UNetwork<>::read("Zakarays_Karate_Club.edges");
	// Create the test data
	auto testData = NetworkManipulator<>::createTestDataRem(refNet, remRatio, seed);
	// Lock test data
	testData.lock();
	auto obsNet = testData.getObsNet();
	// Create an evaluator object
	PerfEvaluator<> perf(testData);
	// Create an ADA predictor
	auto ada = std::make_shared<UADAPredictor<>>(obsNet);
	// Add it to the evaluator
	perf.addPredictor(ada);
	// Create a CNE predictor
	auto cne = std::make_shared<UCNEPredictor<>>(obsNet);
	// Add it to the evaluator
	perf.addPredictor(cne);
	// Create a ROC object
	auto roc = std::make_shared<ROC<>>();
	// Add it to the evaluator
	perf.addPerfMeasure(roc);
	// Create a PR object
	auto pr = std::make_shared<PR<>>();
	// Add it to the evaluator
	perf.addPerfMeasure(pr);
	// Run evaluation
	perf.eval();
	// Print results
	for (auto it = perf.resultsBegin(); it != perf.resultsEnd(); ++it) {
		std::cout << it->first << "\t" << it->second << std::endl;
	}
	return 0;
}
