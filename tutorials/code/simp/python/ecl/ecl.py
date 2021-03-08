# Import the module
import LinkPredPython as lpp
k = 10;
# Create a predictor object
p = lpp.Predictor();
# Load network from file
p.loadnet("Zakarays_Karate_Club.edges");
# Predict top k scores using an encoder-classifier predictor with Node2Vec as encoder and logistic regression as classifier
esv = p.predTopECL(k, "N2V", "LGR");
# Print the scores
print("N2V-LGR");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));
# Predict top k scores using an encoder-classifier predictor with LINE as encoder and feed-forward neural network as classifier
esv = p.predTopECL(k, "LIN", "FFN"); # FFN requires mlpack
# Print the scores
print("LIN-FFN");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));
