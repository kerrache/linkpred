#include <linkpred.hpp>
#include <iostream>

using namespace LinkPred::Simp;

int main() {
	Evaluator eval;
	eval.addADA();
	eval.addPST("PST", "pst.txt"); // Load scores from pst.txt
	eval.addRAL();
	eval.addPR();
	eval.addTPR();
	eval.run("Zakarays_Train.edges", "Zakarays_Test.edges");
	return 0;
}
