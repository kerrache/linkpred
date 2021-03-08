#include <linkpred.hpp>
using namespace LinkPred;
int main() {
	// Read network frome file
	auto net = UNetwork<>::read("Zakarays_Karate_Club.edges");
	// Print to standard output
	net->print();
	return 0;
}
