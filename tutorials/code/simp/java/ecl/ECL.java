public class ECL {
	static {
		// Load the library
		System.loadLibrary("LinkPredJava");
	}
	public static void main(String[] args) {
		int k = 10;
		// Create a prtedictor object
		Predictor p = new Predictor();
		// Load network from file
		p.loadnet("Zakarays_Karate_Club.edges");
		// Predict top k scores using an encoder-classifier predictor with Node2Vec as encoder and logistic regression as classifier
		EdgeScoreVec esv = p.predTopECL(k, "N2V", "LGR");
		// Print scores
		System.out.println("N2V-LGR");
		for (int i = 0; i < esv.size(); i++) {
			EdgeScore es = esv.get(i);
			System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
		} 
		// Predict top k scores using an encoder-classifier predictor with LINE as encoder and feed-forward neural network as classifier
		esv = p.predTopECL(k, "LIN", "FFN"); // FFN requires mlpack
		// Print scores
		System.out.println("LIN-FFN");
		for (int i = 0; i < esv.size(); i++) {
			EdgeScore es = esv.get(i);
			System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
		} 
	}
}
