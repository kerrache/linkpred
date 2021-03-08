# Import the module
import LinkPredPython as lpp
# Create a prtedictor object
p = lpp.Predictor();
# Load network from file
p.loadnet("Zakarays_Karate_Club.edges");
# Compute the score for the two edges (1, 34) and (26,34)
esv = lpp.EdgeScoreVec();
es = lpp.EdgeScore();
es.i = "1";
es.j = "34";
esv.push_back(es);
es.i = "26";
es.j = "34";
esv.push_back(es);
p.predKAB(esv);
# Print the scores
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));
