# Import the module
import LinkPredPython as lpp
nbRuns = 10;
edgeRemRatio = 0.1; 
# Create an evuator object
ev = lpp.Evaluator();
# Add predictors to be evuated
ev.addCNE();
ev.addADA();
ev.addKAB();
# Add performance measures
ev.addROC();
ev.addTPR();
# Run experiment on the specified network
ev.run("Zakarays_Karate_Club.edges", nbRuns, edgeRemRatio);
import sys # For printing
# Print the header row
res = ev.getPerfRes(0);
for r in res:
	sys.stdout.write(r.name + "\t");
sys.stdout.write("\n");
# Print the results of each iteration
for i in range(nbRuns):
	res = ev.getPerfRes(i);
	for r in res:
		sys.stdout.write("{:.4f}".format(r.res) + "\t");
	sys.stdout.write("\n");
