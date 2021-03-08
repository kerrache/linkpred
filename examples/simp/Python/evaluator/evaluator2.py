# Import the library
import LinkPredPython as lpp
import sys # For printing
nbRuns = 10;
# Create an evaluator object
ev = lpp.Evaluator();
ev.addENC(); # Default setting
ev.addENC("ENC-HMSM-LGR", "HMSM", "LGR");
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

