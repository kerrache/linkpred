public class Evaluator5 {
	static {
		System.loadLibrary("LinkPredJava"); // Load the library
	}
	public static void main(String[] args) {
		double edgeRemRatio = 0.1;
		boolean keepConnected = false;
		int seed = 0;
		Evaluator eval = new Evaluator();
		eval.genTestData("Zakarays_Karate_Club.edges", "Zakarays_Train.edges", "Zakarays_Test.edges", edgeRemRatio, keepConnected, seed);
	}
}
