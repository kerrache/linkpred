#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main() {
	// Load test data
	auto testData = NetworkManipulator<>::loadTestDataAdd("net-obs.edges", "net-add.edges");
	std::cout << "Reference network:\n";
	auto refNet = testData.getRefNet();
	refNet->print();
	std::cout << "Observed network:\n";
	auto obsNet = testData.getObsNet();
	obsNet->print();
	std::cout << "Positive examples:" << std::endl;
	for (auto it = testData.posBegin(); it != testData.posEnd(); ++it) {
		auto i = refNet->start(*it);
		auto j = refNet->end(*it);
		std::cout << refNet->getLabel(i) << "\t" << refNet->getLabel(j) << std::endl;
	}
	std::cout << "Negative examples (added edges):" << std::endl;
	for (auto it = testData.negBegin(); it != testData.negEnd(); ++it) {
		auto i = refNet->start(*it);
		auto j = refNet->end(*it);
		std::cout << refNet->getLabel(i) << "\t" << refNet->getLabel(j) << std::endl;
	}
	return 0;
}
