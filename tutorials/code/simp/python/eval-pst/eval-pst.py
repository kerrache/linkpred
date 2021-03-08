# Import the module
import LinkPredPython as lpp
# Create an evuator object
ev = lpp.Evaluator();
# Add predictors to be evuated
ev.addADA();
# Use the method addPST to create a predictor that loads scores from pst.txt
ev.addPST("PST", "pst.txt");
ev.addRAL();
# Add performance measures
ev.addPR();
ev.addTPR();
# Run experiment on the specified network
ev.run("Zakarays_Train.edges", "Zakarays_Test.edges");
