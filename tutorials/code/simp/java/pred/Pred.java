public class Pred {
	static {
		// Load the library
		System.loadLibrary("LinkPredJava");
	}
	public static void main(String[] args) {
		// Create a prtedictor object
		Predictor p = new Predictor();
		// Load network from file
		p.loadnet("Zakarays_Karate_Club.edges");
		// Compute the score for the two edges (1, 34) and (26,34)
		EdgeScoreVec esv = new EdgeScoreVec();
		EdgeScore es;
		es = new EdgeScore();
		es.setI("1");
		es.setJ("34");
		esv.add(es);
		es = new EdgeScore();
		es.setI("26");
		es.setJ("34");
		esv.add(es);
		p.predADA(esv);
		// Print the scores
		for (int i = 0; i < esv.size(); i++) {
			es = esv.get(i);
			System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
		} 
	}
}
