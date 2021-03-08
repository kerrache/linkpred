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
	// Create an instance of the KAB predictor
	auto p = std::make_shared<UKABPredictor<>>(testData.getObsNet());
	// Initialize predictor
	p->init();
	// Train predictor
	p->learn();
	// Create a prediction result object
	auto predResults = std::make_shared<PredResults<>>(testData, p);
	// A map to store prediction results
	PerfResults res;
	// Create a ROC object
	ROC<> roc;
	// Compute ROCAUC
	roc.eval(predResults, res);
	std::cout << "ROCAUC: " << res.at(roc.getName()) << std::endl;
	// Create a PR object
	PR<> pr;
	// Compute PRAUC
	pr.eval(predResults, res);
	std::cout << "PRAUC: " << res.at(pr.getName()) << std::endl;
	// Create a TPR object
	TPR<> tpr(testData.getNbPos());
	// Compute TPR
	tpr.eval(predResults, res);
	std::cout << "TPR: " << res.at(tpr.getName()) << std::endl;
	return 0;
}
