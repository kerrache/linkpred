public class Evaluator1 {
	static {
		System.loadLibrary("LinkPredJava"); // Load the library
	}
	public static void main(String[] args) {
		int nbRuns = 10;
		Evaluator eval = new Evaluator();
		eval.addCNE();
		eval.addADA();
		eval.addKAB("KAB-2", 2); // Horizon limit = 2
		eval.addKAB("KAB-4", 4); // Horizon limit = 4
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
