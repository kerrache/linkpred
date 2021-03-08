# Import the library
import LinkPredPython as lpp
# Create a prtedictor object
p = lpp.Predictor();
# Load network from file
p.loadnet("Zakarays_Karate_Club.edges");
esv = p.predAllKAB();
# Print the scores
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));
