# Import the library
import LinkPredPython as lpp
import sys # For printing
nbRuns = 10;
# Create an evaluator object
ev = lpp.Evaluator();
ev.addADA();
ev.addRAL();
ev.addPR();
ev.addTPR();
ev.run("Zakarays_Karate_Club_Train.edges", "Zakarays_Karate_Club_Test.edges");
res = ev.getPerfRes(0);
for r in res:
	sys.stdout.write(r.name + "\t");
sys.stdout.write("\n");
for r in res:
	sys.stdout.write("{:.4f}".format(r.res) + "\t");
sys.stdout.write("\n");
