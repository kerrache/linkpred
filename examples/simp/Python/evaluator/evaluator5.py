# Import the library
import LinkPredPython as lpp
edgeRemRatio = 0.1;
keepConnected = False;
seed = 0;
# Create an evaluator object
ev = lpp.Evaluator();
ev.genTestData("Zakarays_Karate_Club.edges", "Zakarays_Train.edges", "Zakarays_Test.edges", edgeRemRatio, keepConnected, seed);
