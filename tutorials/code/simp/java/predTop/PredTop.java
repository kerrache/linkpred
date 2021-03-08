public class PredTop {
	static {
		// Load the library
		System.loadLibrary("LinkPredJava");
	}
	public static void main(String [] args) {
		int k = 10; // Find top 10
		// Create a prtedictor object
		Predictor p = new Predictor();
		// Load network from file
		p.loadnet("Zakarays_Karate_Club.edges");
		// Predict the top k edges using Adamic Adar index
		EdgeScoreVec esv = p.predTopADA(k);
		// Print the scores
		for (int i = 0; i < esv.size(); i++) {
			EdgeScore es = esv.get(i);
			System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
		}
	}
}
