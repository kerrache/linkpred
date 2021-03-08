#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;
int main() {
	long int seed = 777;
	// Read network from file
	auto net = UNetwork<>::read("Zakarays_Karate_Club.edges");
	// Create a HMSM encoder
	HMSM<> encoder(net, seed);
	// Set encoding dimension
	encoder.setDim(3);
	// Initialize the encoder
	encoder.init();
	// Embed the network
	encoder.encode();
	// Print node codes
	for (std::size_t i = 0; i < net->getNbNodes(); i++) {
		auto v = encoder.getNodeCode(i);
		std::cout << net->getLabel(i) << "\t";
		for (int k = 0; k < v.size(); k++) {
			std::cout << std::fixed << std::setprecision(4) << v[k] << "\t";
		}
		std::cout << std::endl;
	}
	return 0;
}
