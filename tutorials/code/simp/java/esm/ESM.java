public class ESM {
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
		// Predict top k scores using an encoder-classifier predictor with Node2Vec as encoder and L2 similarity
		EdgeScoreVec esv = p.predTopESM(k, "N2V", "L2");
		// Print scores
		System.out.println("N2V-L2");
		for (int i = 0; i < esv.size(); i++) {
			EdgeScore es = esv.get(i);
			System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
		} 
		// Predict top k scores using an encoder-similarity measure predictor with LINE as encoder and cosine similarity
		esv = p.predTopESM(k, "LIN", "CSM");
		// Print scores
		System.out.println("LIN-CSM");
		for (int i = 0; i < esv.size(); i++) {
			EdgeScore es = esv.get(i);
			System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
		} 
	}
}
