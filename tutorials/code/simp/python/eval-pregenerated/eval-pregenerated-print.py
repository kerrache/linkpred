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
import sys # For printing
res = ev.getPerfRes(0);
for r in res:
	sys.stdout.write(r.name + "\t");
sys.stdout.write("\n");
for r in res:
	sys.stdout.write("{:.4f}".format(r.res) + "\t");
sys.stdout.write("\n");
