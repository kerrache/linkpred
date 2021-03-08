#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main() {
	// Remove 20% of the edges
	double remRatio = 0.2;
	long int seed = 777;
	// Read network from file
	auto net = UNetwork<>::read("net.edges");
	// Create the test data
	auto testData = NetworkManipulator<>::createTestDataRem(net, remRatio, seed);
	std::cout << "Reference network:\n";
	testData.getRefNet()->print();
	std::cout << "Observed network:\n";
	testData.getObsNet()->print();
	std::cout << "Positive examples (removed edges):" << std::endl;
	for (auto it = testData.posBegin(); it != testData.posEnd(); ++it) {
		auto i = net->start(*it);
		auto j = net->end(*it);
		std::cout << net->getLabel(i) << "\t" << net->getLabel(j) << std::endl;
	}
	std::cout << "Negative examples:" << std::endl;
	for (auto it = testData.negBegin(); it != testData.negEnd(); ++it) {
		auto i = net->start(*it);
		auto j = net->end(*it);
		std::cout << net->getLabel(i) << "\t" << net->getLabel(j) << std::endl;
	}
	return 0;
}
