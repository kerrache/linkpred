#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main() {
	// Add 20% of the edges
	double addRatio = 0.2;
	long int seed = 888;
	// Read network from file
	auto net = UNetwork<>::read("net.edges");
	// Create the test data
	auto testData = NetworkManipulator<>::createTestDataAdd(net, addRatio, seed);
	std::cout << "Reference network:\n";
	testData.getRefNet()->print();
	std::cout << "Observed network:\n";
	testData.getObsNet()->print();
	std::cout << "Positive examples:" << std::endl;
	for (auto it = testData.posBegin(); it != testData.posEnd(); ++it) {
		auto i = net->start(*it);
		auto j = net->end(*it);
		std::cout << net->getLabel(i) << "\t" << net->getLabel(j) << std::endl;
	}
	std::cout << "Negative examples (added edges):" << std::endl;
	for (auto it = testData.negBegin(); it != testData.negEnd(); ++it) {
		auto i = net->start(*it);
		auto j = net->end(*it);
		std::cout << net->getLabel(i) << "\t" << net->getLabel(j) << std::endl;
	}
	return 0;
}
