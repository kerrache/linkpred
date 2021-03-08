public class Evaluator2 {
	static {
		System.loadLibrary("LinkPredJava"); // Load the library
	}
	public static void main(String[] args) {
		int nbRuns = 10;
		Evaluator eval = new Evaluator();
		eval.addENC(); // Default setting
		eval.addENC("ENC-HMSM-LGR", "HMSM", "LGR");
		eval.addROC();
		eval.addTPR();
		eval.run("Zakarays_Karate_Club.edges", nbRuns);
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
