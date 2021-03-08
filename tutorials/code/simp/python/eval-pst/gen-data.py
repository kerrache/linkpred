# Import the module
import LinkPredPython as lpp
# We remove 10% of the edges
edgeRemRatio = 0.1;
# We will not keep the network connected when removing edges
keepConnected = False;
# Seed of the random number generator
seed = 0;
# Create an Evaluator object
ev = lpp.Evaluator();
# The ground truth network "Zakarays_Karate_Club.edges" is split into an observed network stored in "Zakarays_Train.edges" and a list of removed edges stored in "Zakarays_Test.edges"
ev.genTestData("Zakarays_Karate_Club.edges", "Zakarays_Train.edges", "Zakarays_Test.edges", edgeRemRatio, keepConnected, seed);
