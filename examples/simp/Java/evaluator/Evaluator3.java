public class Evaluator3 {
	static {
		System.loadLibrary("LinkPredJava"); // Load the library
	}
	public static void main(String[] args) {
		int nbRuns = 10;
		double edgeRemRatio = 0.1; 
		// Create an evaluator object
		Evaluator eval = new Evaluator();
		// Add predictors to be evaluated
		eval.addCNE();
		eval.addADA();
		eval.addKAB();
		// Add performance measures
		eval.addROC();
		eval.addTPR();
		// Run experiment on the specified network
		eval.run("Zakarays_Karate_Club.edges", nbRuns, edgeRemRatio);
	}
}
