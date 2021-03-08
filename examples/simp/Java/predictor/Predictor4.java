public class Predictor4 {
	static {
		System.loadLibrary("LinkPredJava"); // Load the library
	}
	public static void main(String[] args) {
		Predictor p = new Predictor(); // Create a predictor object
		p.loadnet("Zakarays_Karate_Club.edges"); // Load network from file
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
		p.predKAB(esv);
		for (int i = 0; i < esv.size(); i++) {
			es = esv.get(i);
			System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
		} 
	}
}
