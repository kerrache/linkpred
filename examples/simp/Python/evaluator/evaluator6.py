# Import the library
import LinkPredPython as lpp
# Create an evuator object
ev = lpp.Evaluator();
ev.addADA();
# Load scores from pst.txt
ev.addPST("PST", "pst.txt");
ev.addRAL();
ev.addPR();
ev.addTPR();
ev.run("Zakarays_Train.edges", "Zakarays_Test.edges");
