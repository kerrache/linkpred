public class GenData {
	static {
		// Load the library
		System.loadLibrary("LinkPredJava");
	}
	public static void main(String[] args) {
		// We remove 10% of the edges
		double edgeRemRatio = 0.1;
		// We will not keep the network connected when removing edges
		boolean keepConnected = false;
		// Seed of the random number generator
		int seed = 0;
		// Create an Evaluator object
		Evaluator eval = new Evaluator();
		// The ground truth network "Zakarays_Karate_Club.edges" is split into an observed network stored in "Zakarays_Train.edges" and a list of removed edges stored in "Zakarays_Test.edges"
		eval.genTestData("Zakarays_Karate_Club.edges", "Zakarays_Train.edges", "Zakarays_Test.edges", edgeRemRatio, keepConnected, seed);
	}
}
