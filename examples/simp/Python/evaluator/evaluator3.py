# Import the library
import LinkPredPython as lpp
nbRuns = 10;
edgeRemRatio = 0.1;
# Create an evaluator object
ev = lpp.Evaluator();
# Add predictors to be evaluated
ev.addCNE();
ev.addADA();
ev.addKAB();
# Add performance measures
ev.addROC();
ev.addTPR();
# Run experiment on the specified network
ev.run("Zakarays_Karate_Club.edges", nbRuns, edgeRemRatio);
