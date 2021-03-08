class MEdgeScore:
	i = "";
	j = "";
	score = 0;

# Import the library
import LinkPredPython as lpp
# Create a predictor object
p = lpp.Predictor();
# Load network from file
p.loadnet("Zakarays_Karate_Club.edges");
# Predict the score of all non-exisitng edges using Adamic Adar index
esv = p.predAllADA();
mes = MEdgeScore();
# Print the scores
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));
	mes.i = es.i;
	mes.j = es.j;
	mes.score = es.score;
	print(mes.i + "\t" + mes.j + "\t" + "{:.4f}".format(mes.score));
	
	
