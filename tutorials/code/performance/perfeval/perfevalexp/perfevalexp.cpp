#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
// This is a factory to create predictors and performance measures
class Factory: public PEFactory<> {
public:
	// This method creates predictors
	virtual std::vector<std::shared_ptr<ULPredictor<>>> getPredictors(std::shared_ptr<UNetwork<> const> obsNet) {
		std::vector<std::shared_ptr<ULPredictor<>>> prs;
		// Add ADA
		prs.push_back(std::make_shared<URALPredictor<>>(obsNet));
		// Add KAB
		prs.push_back(std::make_shared<UKABPredictor<>>(obsNet));
		return prs;
	}
	// This method creates performance measures
	virtual std::vector<std::shared_ptr<PerfMeasure<>>> getPerfMeasures(TestData<> const & testData) {
		std::vector<std::shared_ptr<PerfMeasure<>>> pms;
		// Add PR
		pms.push_back(std::make_shared<PR<>>());
		// Add ROC
		pms.push_back(std::make_shared<ROC<>>());
		return pms;
	}
	virtual ~Factory() = default;
};
int main() {
	// Remove 10% of the edges
	double remRatio = 0.1;
	long int seed = 888;
	// Read network from file
	auto refNet = UNetwork<>::read("Zakarays_Karate_Club.edges");
	// The parameters of our experiment
	PerfeEvalExpDescp<> ped;
	// Set the reference network
	ped.refNet = refNet;
	// We run 10 tests
	ped.nbTestRuns = 10;
	ped.seed = 777;
	// Create the factory
	auto factory = std::make_shared<Factory>();
	// Create the experiment object
	PerfEvalExp<> exp(ped, factory);
	// Run the experiment 
	exp.run();
	return 0;
}
