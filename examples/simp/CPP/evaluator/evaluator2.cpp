#include <linkpred.hpp>
#include <iostream>

using namespace LinkPred::Simp;

int main() {
	int nbRuns = 10;
	Evaluator eval;
	eval.addECL(); // Default setting
	eval.addECL("ENC-HMSM-LGR", "HMSM", "LGR");
	eval.addROC();
	eval.addTPR();
	eval.run("Zakarays_Karate_Club.edges", nbRuns);
	for (int i = 0; i < nbRuns; i++) {
		auto res = eval.getPerfRes(i);
		for (auto it = res.begin(); it != res.end(); ++it) {
			std::cout << it->name << "\t" << it->res << "\t";
		}
		std::cout << "\n";
	}
	return 0;
}
