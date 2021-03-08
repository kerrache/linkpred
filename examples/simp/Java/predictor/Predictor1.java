public class Predictor1 {
	static {
		System.loadLibrary("LinkPredJava"); // Load the library
	}
	public static void main(String[] args) {
		Predictor p = new Predictor(); // Create a predictor object
		p.loadnet("Zakarays_Karate_Club.edges"); // Load network from file
		EdgeScoreVec esv = p.predAllADA(); // Predict the score of all non-exisitng edges using Adamic Adar index.
		for (int i = 0; i < esv.size(); i++) {
			EdgeScore es = esv.get(i);
			System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
		} 
	}
}
