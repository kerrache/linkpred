public class EvalAutoPrint {
	static {
		// Load the library
		System.loadLibrary("LinkPredJava");
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
		// Print the header row
		PerfResVec res = eval.getPerfRes(0);
		for (int j = 0; j < res.size(); j++) {
			System.out.print(res.get(j).getName() + "\t") ;
		}
		System.out.println();
		// Print the results of each iteration
		for(int i = 0; i < nbRuns; i++) {
			res = eval.getPerfRes(i);
			for (int j = 0; j < res.size(); j++) {
				System.out.printf("%.4f\t", res.get(j).getRes()) ;
			}
			System.out.println();
		}
	}
}
