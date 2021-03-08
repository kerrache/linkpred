# Import the library
import LinkPredPython as lpp
import sys # For printing
nbRuns = 10;
# Create an evaluator object
ev = lpp.Evaluator();
ev.addCNE();
ev.addADA();
ev.addKAB("KAB-2", 2); # Horizon limit = 2
ev.addKAB("KAB-4", 4); # Horizon limit = 4
ev.addROC();
ev.addTPR();
ev.run("Zakarays_Karate_Club.edges", nbRuns);
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

