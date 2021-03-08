public class Evaluator4 {
	static {
		System.loadLibrary("LinkPredJava"); // Load the library
	}
	public static void main(String[] args) {
		Evaluator eval = new Evaluator();
		eval.addADA();
		eval.addRAL();
		eval.addPR();
		eval.addTPR();
		eval.run("Zakarays_Karate_Club_Train.edges", "Zakarays_Karate_Club_Test.edges");
		PerfResVec res = eval.getPerfRes(0);
		for (int j = 0; j < res.size(); j++) {
			System.out.print(res.get(j).getName() + "\t") ;
		}
		System.out.println();
		for (int j = 0; j < res.size(); j++) {
			System.out.printf("%.4f\t", res.get(j).getRes()) ;
		}
		System.out.println();
	}
}
