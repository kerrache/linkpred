import LinkPredPython as lpp
k = 10;
p = lpp.Predictor();
p.loadnet("Zakarays_Karate_Club.edges");

esv = p.predTopENC(k, "DPW", "FFN");
print("Top " + str(k) + " using DPW-FFN");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "DPW", "LSVM");
print("Top " + str(k) + " using DPW-LSVM");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "DPW", "LGR");
print("Top " + str(k) + " using DPW-LGR");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "DPW", "NVB");
print("Top " + str(k) + " using DPW-NVB");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));


esv = p.predTopENC(k, "HMSM", "FFN");
print("Top " + str(k) + " using HMSM-FFN");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "HMSM", "LSVM");
print("Top " + str(k) + " using HMSM-LSVM");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "HMSM", "LGR");
print("Top " + str(k) + " using HMSM-LGR");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "HMSM", "NVB");
print("Top " + str(k) + " using HMSM-NVB");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));


esv = p.predTopENC(k, "LVS", "FFN");
print("Top " + str(k) + " using LVS-FFN");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "LVS", "LSVM");
print("Top " + str(k) + " using LVS-LSVM");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "LVS", "LGR");
print("Top " + str(k) + " using LVS-LGR");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "LVS", "NVB");
print("Top " + str(k) + " using LVS-NVB");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));


esv = p.predTopENC(k, "LIN", "FFN");
print("Top " + str(k) + " using LIN-FFN");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "LIN", "LSVM");
print("Top " + str(k) + " using LIN-LSVM");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "LIN", "LGR");
print("Top " + str(k) + " using LIN-LGR");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "LIN", "NVB");
print("Top " + str(k) + " using LIN-NVB");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));


esv = p.predTopENC(k, "N2V", "FFN");
print("Top " + str(k) + " using N2V-FFN");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "N2V", "LSVM");
print("Top " + str(k) + " using N2V-LSVM");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "N2V", "LGR");
print("Top " + str(k) + " using N2V-LGR");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

esv = p.predTopENC(k, "N2V", "NVB");
print("Top " + str(k) + " using N2V-NVB");
for es in esv:
	print(es.i + "\t" + es.j + "\t" + "{:.4f}".format(es.score));

