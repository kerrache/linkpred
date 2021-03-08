public class Predictor3 {
	static {
		System.loadLibrary("LinkPredJava"); // Load the library
	}
	public static void main(String[] args) {
		int k = 10;
		Predictor p = new Predictor(); // Create a predictor object
		p.loadnet("Zakarays_Karate_Club.edges"); // Load network from file
		{
			EdgeScoreVec esv = p.predTopECL(k, "DPW", "FFN");
			System.out.println("Top " + k + " using DPW-FFN");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "DPW", "LSVM");
			System.out.println("Top " + k + " using DPW-LSVM");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "DPW", "LGR");
			System.out.println("Top " + k + " using DPW-LGR");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "DPW", "NVB");
			System.out.println("Top " + k + " using DPW-NVB");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}

		{
			EdgeScoreVec esv = p.predTopECL(k, "HMSM", "FFN");
			System.out.println("Top " + k + " using HMSM-FFN");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "HMSM", "LSVM");
			System.out.println("Top " + k + " using HMSM-LSVM");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "HMSM", "LGR");
			System.out.println("Top " + k + " using HMSM-LGR");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "HMSM", "NVB");
			System.out.println("Top " + k + " using HMSM-NVB");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}

		{
			EdgeScoreVec esv = p.predTopECL(k, "LVS", "FFN");
			System.out.println("Top " + k + " using LVS-FFN");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "LVS", "LSVM");
			System.out.println("Top " + k + " using LVS-LSVM");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "LVS", "LGR");
			System.out.println("Top " + k + " using LVS-LGR");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "LVS", "NVB");
			System.out.println("Top " + k + " using LVS-NVB");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}

		{
			EdgeScoreVec esv = p.predTopECL(k, "LIN", "FFN");
			System.out.println("Top " + k + " using LIN-FFN");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "LIN", "LSVM");
			System.out.println("Top " + k + " using LIN-LSVM");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "LIN", "LGR");
			System.out.println("Top " + k + " using LIN-LGR");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "LIN", "NVB");
			System.out.println("Top " + k + " using LIN-NVB");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}

		{
			EdgeScoreVec esv = p.predTopECL(k, "N2V", "FFN");
			System.out.println("Top " + k + " using N2V-FFN");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "N2V", "LSVM");
			System.out.println("Top " + k + " using N2V-LSVM");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "N2V", "LGR");
			System.out.println("Top " + k + " using N2V-LGR");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
		{
			EdgeScoreVec esv = p.predTopECL(k, "N2V", "NVB");
			System.out.println("Top " + k + " using N2V-NVB");
			for (int i = 0; i < esv.size(); i++) {
				EdgeScore es = esv.get(i);
				System.out.println(es.getI() + "\t" + es.getJ() + "\t" + es.getScore());
			}
		}
	}
}
