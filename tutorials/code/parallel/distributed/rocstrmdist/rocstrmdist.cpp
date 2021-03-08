#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main(int argc, char*argv[]) {
	double remRatio = 0.1;
	long int seed = 777;
	// Initialize MPI
	MPI_Init(&argc, &argv);
	int procID = 0;
	// Get local process ID
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);
	// Read network from file
	auto net = UNetwork<>::read("Infectious.edges");
	// Create the test data
	auto testData = NetworkManipulator<>::createTestDataRem(net, remRatio, seed, false);
	testData.lock();
	// Create an ADA predictor
	auto p = std::make_shared<UADAPredictor<>>(testData.getObsNet());
	// Initialize predictor
	p->init();
	// Train predictor
	p->learn();
	// Create object to store prediction results
	auto predResults = std::make_shared<PredResults<>>(testData, p);
	// Create a ROC object
	auto roc = std::make_shared<ROC<>>("ROC");
	// Set options
	roc->setParallel(true);
	roc->setDistributed(true);
	roc->setStrmEnabled(true);
	// Create object to store results
	PerfResults res;
	// Evaluate performance
	roc->eval(predResults, res);
	// Print results
	if (procID == 0) {
		std::cout << "#ROCAUC (streaming): " << res.at(roc->getName()) << std::endl;
	}
	// Finalize MPI
	MPI_Finalize();
	return 0;
}
