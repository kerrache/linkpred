public class EvalPregeneratedPrint {
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
		// Print the header row
		PerfResVec res = eval.getPerfRes(0);
		for (int j = 0; j < res.size(); j++) {
			System.out.print(res.get(j).getName() + "\t") ;
		}
		System.out.println();
		// Print the results
		for (int j = 0; j < res.size(); j++) {
			System.out.printf("%.4f\t", res.get(j).getRes()) ;
		}
		System.out.println();
	}
}
