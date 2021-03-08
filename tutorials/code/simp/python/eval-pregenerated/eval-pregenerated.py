# Import the module
import LinkPredPython as lpp
# Create an evuator object
ev = lpp.Evaluator();
# Add predictors to be evuated
ev.addADA();
ev.addRAL();
# Add performance measures
ev.addPR();
ev.addTPR();
# Run experiment on the specified network
ev.run("Zakarays_Karate_Club_Train.edges", "Zakarays_Karate_Club_Test.edges");
