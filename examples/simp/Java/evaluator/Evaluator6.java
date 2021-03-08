public class Evaluator6 {
	static {
		System.loadLibrary("LinkPredJava"); // Load the library
	}
	public static void main(String[] args) {
		Evaluator eval = new Evaluator();
		eval.addADA();
		// Load scores from pst.txt
		eval.addPST("PST", "pst.txt");
		eval.addRAL();
		eval.addPR();
		eval.addTPR();
		eval.run("Zakarays_Train.edges", "Zakarays_Test.edges");
	}
}
