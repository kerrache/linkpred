#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main() {
	auto net = UNetwork<>::read("Infectious.edges");
	UCNEPredictor<> predictor(net);
	predictor.init();
	predictor.learn();
	std::cout << "#Start\tEnd\tScore\n";
	for (auto it=net->nonEdgesBegin();it!=net->nonEdgesEnd();++it){
		auto i = net->getLabel(net->start(*it));
		auto j = net->getLabel(net->end(*it));
		double sc = predictor.score(*it);
		std::cout << i << "\t" << j << "\t" << sc << std::endl;
	}
	return 0;
}
