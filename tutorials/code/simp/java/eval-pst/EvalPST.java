public class EvalPST {
	static {
		// Load the library
		System.loadLibrary("LinkPredJava");
	}
	public static void main(String[] args) {
		// Create an evaluator object
		Evaluator eval = new Evaluator();
		// Add predictors to be evaluated
		eval.addADA();
		// Use the method addPST to create a predictor that loads scores from pst.txt
		eval.addPST("PST", "pst.txt");
		eval.addRAL();
		// Add performance measures
		eval.addPR();
		eval.addTPR();
		// Run experiment on the specified network
		eval.run("Zakarays_Train.edges", "Zakarays_Test.edges");
	}
}
