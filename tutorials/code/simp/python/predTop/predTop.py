# Import the module
import LinkPredPython as lpp
k = 10;
# Create a predictor object
p = lpp.Predictor();
# Load network from file
p.loadnet("Zakarays_Karate_Club.edges");
# Predict the top k edges using Adamic Adar index
esv = p.predTopADA(k); 
# Print the scores
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));
