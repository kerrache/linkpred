public class Main {
	static {
		System.loadLibrary("LinkPredJava"); // Load the required library
	}
	public static void main(String[] args) {
		/*
		Predictor p = new Predictor(); // Create a predictor object
		p.loadnet("karate.edges"); // Read the network
		EdgeScoreVec esv = p.cneAll(); // Predict the score for all negative edges
		for (int i = 0; i < esv.size(); i++) {
			EdgeScore es = esv.get(i);
			System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
		} 
		*/
		int nbRuns = 10;
		Evaluator eval = new Evaluator(); // Create an evaluator object
		eval.addCNE();
		eval.addCRA();
		eval.addROC();
		eval.addTPR();
		eval.run("Zakarays_Karate_Club.edges", nbRuns);
		for(int i = 0; i < nbRuns; i++) {
			PerfResVec res = eval.getPerfRes(i);
			for (int j = 0; j < res.size(); j++) {
				System.out.print(res.get(j).getName() + "\t");
				System.out.printf("%.4f", res.get(j).getRes());
				System.out.print("\t");
			}
			System.out.println();
		} 

	}
}
