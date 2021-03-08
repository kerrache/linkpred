public class EvalPregenerated {
	static {
		// Load the library
		System.loadLibrary("LinkPredJava");
	}
	public static void main(String[] args) {
		// Create an evaluator object
		Evaluator eval = new Evaluator();
		// Add predictors to be evaluated
		eval.addADA();
		eval.addRAL();
		// Add performance measures
		eval.addPR();
		eval.addTPR();
		// Run experiment on the specified train and test set
		eval.run("Zakarays_Karate_Club_Train.edges", "Zakarays_Karate_Club_Test.edges");
	}
}
