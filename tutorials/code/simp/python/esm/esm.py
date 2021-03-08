# Import the module
import LinkPredPython as lpp
k = 10;
# Create a predictor object
p = lpp.Predictor();
# Load network from file
p.loadnet("Zakarays_Karate_Club.edges");
# Predict top k scores using an encoder-classifier predictor with Node2Vec as encoder and L2 similarity
esv = p.predTopESM(k, "N2V", "L2");
# Print the scores
print("N2V-L2");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));
# Predict top k scores using an encoder-similarity measure predictor with LINE as encoder and cosine similarity
esv = p.predTopESM(k, "LIN", "CSM");
# Print the scores
print("LIN-CSM");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));
