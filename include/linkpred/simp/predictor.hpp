/*
 * This file is part of LinkPred.
 *
 * LinkPred: A high performance library for link prediction in complex networks.
 * Copyright (C) 2021  by Said Kerrache.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * \file
 * @ingroup Simp
 * @brief Contains the definition of a class that simplifies the use of link prediction algorithms.
 */

#ifndef SIMPPREDICTOR_HPP_
#define SIMPPREDICTOR_HPP_

#include "LinkPredConfig.hpp"
#include "linkpred/simp/edgescore.hpp"
#include "linkpred/predictors/predictors.hpp"
#include "linkpred/graphalg/encoders/encoders.hpp"
#include "linkpred/ml/classifiers/classifiers.hpp"
#include "linkpred/ml/simmeasures/simmeasures.hpp"

namespace LinkPred {
namespace Simp {

/**
 * @brief A class that simplifies the use of link prediction algorithms.
 */
class Predictor {
private:
	std::shared_ptr<UNetwork<>> net; /**< The network.*/

#ifdef LINKPRED_WITH_OPENMP
	bool parallel = false; /**< Enable/disable parallelism. */
#endif

	/**
	 * @param predictor A link prediction algorithm.
	 * @return The score of all negative links computed using predictor.
	 */
	std::vector<EdgeScore> predAll(std::shared_ptr<ULPredictor<>> predictor);

	/**
	 * @param predictor A link prediction algorithm. Node IDs are used instead of labels.
	 * @return The score of all negative links computed using predictor.
	 */
	std::vector<EdgeScoreByID> predAllByID(
			std::shared_ptr<ULPredictor<>> predictor);

	/**
	 * @param predictor A link prediction algorithm.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 */
	void pred(std::shared_ptr<ULPredictor<>> predictor,
			std::vector<EdgeScore> &edgeScores);

	/**
	 * @param predictor A link prediction algorithm.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTop(std::shared_ptr<ULPredictor<>> predictor,
			int k);

	/**
	 * Create an encoder based on its name.
	 * @param name The name of the encoder.
	 * @param net The network.
	 * @param dim The dimension of the embedding.
	 * @param seed Seed of the random number generator.
	 * @return The encoder.
	 */
	std::shared_ptr<Encoder<>> createEncoder(std::string name,
			std::shared_ptr<const UNetwork<>> net, int dim, long int seed);

	/**
	 * Create a classifier based on its name.
	 * @param name The name of the classifier.
	 * @param dim The dimension of the input (edge embedding).
	 * @param seed Seed of the random number generator.
	 * @return The classifier.
	 */
	std::shared_ptr<Classifier<>> createClassifier(std::string name, int dim,
			long int seed);

	/**
	 * Create a similarity measure based on its name.
	 * @param name The name of the similarity measure.
	 * @param seed Seed of the random number generator.
	 * @return The similarity measure.
	 */
	std::shared_ptr<SimMeasure> createSimMeasure(std::string name);

public:

	/**
	 * @param net The network.
	 */
	Predictor() = default;

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	Predictor(Predictor const &that) = default;

	/**
	 * Destructor.
	 */
	virtual ~Predictor() = default;

	/**
	 * Returns the number of nodes in the network.
	 * @return The number of nodes in the network.
	 */
	int getNbNodes() const;

	/**
	 * Translate node ID to label.
	 * @param i Node internal ID (sequential from 0 to nbNodes-1).
	 * @return The node label.
	 */
	std::string getLabel(int i) const;

	/**
	 * Translate node label to ID.
	 * @param i Node label.
	 * @return The node ID.
	 */
	int getID(std::string const &i) const;

	/**
	 * Check if an edge exists using internal node IDs.
	 * @param i ID of start node.
	 * @param j ID of end node.
	 * @return True if (i, j) is an edge, false otherwise.
	 */
	bool isEdgeByID(int i, int j) const;

	/**
	 * Check if an edge exists using internal node IDs.
	 * @param i ID of start node.
	 * @param j ID of end node.
	 * @return True if (i, j) is an edge, false otherwise.
	 */
	bool isEdgeByLabel(std::string const &i, std::string const &j) const;

	/**
	 * Load network from file. The format is list of edges. Comments must be on separate lines and start with #.
	 * @param fileName The file name.
	 */
	void loadnet(std::string fileName);

	/**
	 * Adamic Adar predictor.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllADA();

	/**
	 * Adamic Adar predictor.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 */
	void predADA(std::vector<EdgeScore> &edgeScores);

	/**
	 * Adamic Adar predictor.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopADA(int k);

	/**
	 * Common neighbors.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllCNE();

	/**
	 * Common neighbors.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 */
	void predCNE(std::vector<EdgeScore> &edgeScores);

	/**
	 * Common neighbors.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopCNE(int k);

	/**
	 * Cannistraci resource allocation.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllCRA();

	/**
	 * Cannistraci resource allocation.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 */
	void predCRA(std::vector<EdgeScore> &edgeScores);

	/**
	 * Cannistraci resource allocation.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopCRA(int k);

	/**
	 * Encoder-classifier link predictor.
	 * @param encoderName The name of the encoder. Possible values are: DPW (DeepWalk), HMSM (Hidden Metric Space Model), LVS (LargeVis), LEM (Laplacian Eigenmaps), LIN (LINE), LLE (Locally Linear Embedding), MFC (Matrix Factorization), and N2V (Node2Vec).
	 * @param classifierName The name of the classifier. Possible values are: FFN (feed-forward neural network withn default architecture), LSVM (linear SVM), LGR (logistic regression), NVB (naive Bayes). All classifiers except logistic regression requirte compilation with mlpack.
	 * @param dim The dimension of the embedding space. If set to zero, the default value is used (the default dimension depends on the ecnoder).
	 * @param posRatio Ratio of positive edges used in the training of the classifier.
	 * @param negRatio Ratio of negative edges used in the training of the classifier.
	 * @param seed Seed of the random number generator.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllECL(std::string encoderName = "N2V",
			std::string classifierName = "LGR", int dim = 0, double posRatio =
					1.0, double negRatio = 1.0, long int seed = 0);

	/**
	 * Encoder-classifier link predictor.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 * @param encoderName The name of the encoder. Possible values are: DPW (DeepWalk), HMSM (Hidden Metric Space Model), LVS (LargeVis), LEM (Laplacian Eigenmaps), LIN (LINE), LLE (Locally Linear Embedding), MFC (Matrix Factorization), and N2V (Node2Vec).
	 * @param classifierName The name of the classifier. Possible values are: FFN (feed-forward neural network withn default architecture), LSVM (linear SVM), LGR (logistic regression), NVB (naive Bayes). All classifiers except logistic regression requirte compilation with mlpack.
	 * @param dim The dimension of the embedding space. If set to zero, the default value is used (the default dimension depends on the ecnoder).
	 * @param posRatio Ratio of positive edges used in the training of the classifier.
	 * @param negRatio Ratio of negative edges used in the training of the classifier.
	 * @param seed Seed of the random number generator.
	 */
	void predECL(std::vector<EdgeScore> &edgeScores, std::string encoderName =
			"N2V", std::string classifierName = "LGR", int dim = 0,
			double posRatio = 1.0, double negRatio = 1.0, long int seed = 0);

	/**
	 * Encoder-classifier link predictor.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @param encoderName The name of the encoder. Possible values are: DPW (DeepWalk), HMSM (Hidden Metric Space Model), LVS (LargeVis),LEM (Laplacian Eigenmaps), LIN (LINE), LLE (Locally Linear Embedding), MFC (Matrix Factorization), and N2V (Node2Vec).
	 * @param classifierName The name of the classifier. Possible values are: FFN (feed-forward neural network withn default architecture), LSVM (linear SVM), LGR (logistic regression), NVB (naive Bayes). All classifiers except logistic regression requirte compilation with mlpack.
	 * @param dim The dimension of the embedding space. If set to zero, the default value is used (the default dimension depends on the ecnoder).
	 * @param posRatio Ratio of positive edges used in the training of the classifier.
	 * @param negRatio Ratio of negative edges used in the training of the classifier.
	 * @param seed Seed of the random number generator.
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopECL(int k, std::string encoderName = "N2V",
			std::string classifierName = "LGR", int dim = 0, double posRatio =
					1.0, double negRatio = 1.0, long int seed = 0);

	/**
	 * Encoder-similarity measure link predictor.
	 * @param encoderName The name of the encoder. Possible values are: DPW (DeepWalk), HMSM (Hidden Metric Space Model), LVS (LargeVis), LEM (Laplacian Eigenmaps), LIN (LINE), LLE (Locally Linear Embedding), MFC (Matrix Factorization), and N2V (Node2Vec).
	 * @param simMeasureName The name of the similarity measure. Possible values are: CSM (cosine similarity), DTP (dot product), L1 (L1 similarity), L2 (L2 similarity), PRS (Pearson similarity).
	 * @param dim The dimension of the embedding space. If set to zero, the default value is used (the default dimension depends on the ecnoder).
	 * @param seed Seed of the random number generator.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllESM(std::string encoderName = "N2V",
			std::string simMeasureName = "L2", int dim = 0, long int seed = 0);

	/**
	 * Encoder-similarity measure link predictor.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 * @param encoderName The name of the encoder. Possible values are: DPW (DeepWalk), HMSM (Hidden Metric Space Model), LVS (LargeVis),LEM (Laplacian Eigenmaps), LIN (LINE), LLE (Locally Linear Embedding), MFC (Matrix Factorization), and N2V (Node2Vec).
	 * @param simMeasureName The name of the similarity measure. Possible values are: CSM (cosine similarity), DTP (dot product), L1 (L1 similarity), L2 (L2 similarity), PRS (Pearson similarity).
	 * @param dim The dimension of the embedding space. If set to zero, the default value is used (the default dimension depends on the ecnoder).
	 * @param seed Seed of the random number generator.
	 */
	void predESM(std::vector<EdgeScore> &edgeScores, std::string encoderName =
			"N2V", std::string simMeasureName = "L2", int dim = 0,
			long int seed = 0);

	/**
	 * Encoder-similarity measure link predictor.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @param encoderName The name of the encoder. Possible values are: DPW (DeepWalk), HMSM (Hidden Metric Space Model), LVS (LargeVis), LEM (Laplacian Eigenmaps), LIN (LINE), LLE (Locally Linear Embedding), MFC (Matrix Factorization), and N2V (Node2Vec).
	 * @param simMeasureName The name of the similarity measure. Possible values are: CSM (cosine similarity), DTP (dot product), L1 (L1 similarity), L2 (L2 similarity), PRS (Pearson similarity).
	 * @param dim The dimension of the embedding space. If set to zero, the default value is used (the default dimension depends on the ecnoder).
	 * @param seed Seed of the random number generator.
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopESM(int k, std::string encoderName = "N2V",
			std::string simMeasureName = "L2", int dim = 0, long int seed = 0);

	/**
	 * Fast blocking model.
	 * @param maxIter Maximum number of iterations.
	 * @param seed Seed of the random number generator.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllFBM(int maxIter = 50, long int seed = 0);

	/**
	 * Fast blocking model.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 * @param maxIter Maximum number of iterations.
	 * @param seed Seed of the random number generator.
	 */
	void predFBM(std::vector<EdgeScore> &edgeScores, int maxIter = 50,
			long int seed = 0);

	/**
	 * Fast blocking model.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @param maxIter Maximum number of iterations.
	 * @param seed Seed of the random number generator.
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopFBM(int k, int maxIter = 50,
			long int seed = 0);

	/**
	 * Hub depromoted index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllHDI();

	/**
	 * Hub depromoted index.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 */
	void predHDI(std::vector<EdgeScore> &edgeScores);

	/**
	 * Hub depromoted index.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopHDI(int k);

	/**
	 * Hub promoted index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllHPI();

	/**
	 * Hub promoted index.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 */
	void predHPI(std::vector<EdgeScore> &edgeScores);

	/**
	 * Hub promoted index.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopHPI(int k);

	/**
	 * Hierarchical random graph.
	 * @param nbBeans Number of bins in edge statistics histogram.
	 * @param nbSamples Number of samples to take for predictions.
	 * @param seed Seed of the random number generator.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllHRG(int nbBeans = 25, int nbSamples = 10000,
			long int seed = 0);

	/**
	 * Hierarchical random graph.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 * @param nbBeans Number of bins in edge statistics histogram.
	 * @param nbSamples Number of samples to take for predictions.
	 * @param seed Seed of the random number generator.
	 */
	void predHRG(std::vector<EdgeScore> &edgeScores, int nbBeans = 25,
			int nbSamples = 10000, long int seed = 0);

	/**
	 * Hierarchical random graph.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @param nbBeans Number of bins in edge statistics histogram.
	 * @param nbSamples Number of samples to take for predictions.
	 * @param seed Seed of the random number generator.
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopHRG(int k, int nbBeans = 25, int nbSamples =
			10000, long int seed = 0);

	/**
	 * Hypermap.
	 * @param m The parameter m (see the algorithm description).
	 * @param L The parameter L (see the algorithm description).
	 * @param gamma The power law exponent gamma (see the algorithm description).
	 * @param zeta The parameter zeta (see the algorithm description).
	 * @param T The parameter T (see the algorithm description).
	 * @param seed The random number generator seed.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllHYP(double m = 1.5, double L = 1,
			double gamma = 2.1, double zeta = 1, double T = 0.8, long int seed =
					0);

	/**
	 * Hypermap.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 * @param m The parameter m (see the algorithm description).
	 * @param L The parameter L (see the algorithm description).
	 * @param gamma The power law exponent gamma (see the algorithm description).
	 * @param zeta The parameter zeta (see the algorithm description).
	 * @param T The parameter T (see the algorithm description).
	 * @param seed The random number generator seed.
	 */
	void predHYP(std::vector<EdgeScore> &edgeScores, double m = 1.5, double L =
			1, double gamma = 2.1, double zeta = 1, double T = 0.8,
			long int seed = 0);

	/**
	 * Hypermap.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @param m The parameter m (see the algorithm description).
	 * @param L The parameter L (see the algorithm description).
	 * @param gamma The power law exponent gamma (see the algorithm description).
	 * @param zeta The parameter zeta (see the algorithm description).
	 * @param T The parameter T (see the algorithm description).
	 * @param seed The random number generator seed.
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopHYP(int k, double m = 1.5, double L = 1,
			double gamma = 2.1, double zeta = 1, double T = 0.8, long int seed =
					0);

	/**
	 * Jackard index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllJID();

	/**
	 * Jackard index.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 */
	void predJID(std::vector<EdgeScore> &edgeScores);

	/**
	 * Jackard index.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopJID(int k);

	/**
	 * A scalable popularity similarity link predictor proposed in: "Kerrache, S., Alharbi,
	 * R. & Benhidour, H. A Scalable Similarity-Popularity Link Prediction Method. Sci Rep 10, 6394 (2020)".
	 * @param horizLim Horizon limit.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllKAB(int horizLim = 2);

	/**
	 * A scalable popularity similarity link predictor proposed in: "Kerrache, S., Alharbi,
	 * R. & Benhidour, H. A Scalable Similarity-Popularity Link Prediction Method. Sci Rep 10, 6394 (2020)".
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 * @param horizLim Horizon limit.
	 */
	void predKAB(std::vector<EdgeScore> &edgeScores, int horizLim = 2);

	/**
	 * A scalable popularity similarity link predictor proposed in: "Kerrache, S., Alharbi,
	 * R. & Benhidour, H. A Scalable Similarity-Popularity Link Prediction Method. Sci Rep 10, 6394 (2020)".
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @param horizLim Horizon limit.
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopKAB(int k, int horizLim = 2);

	/**
	 * Local path.
	 * @param epsilon The weight of paths of length 3.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllLCP(double epsilon = 0.001);

	/**
	 * Local path.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 * @param epsilon The weight of paths of length 3.
	 */
	void predLCP(std::vector<EdgeScore> &edgeScores, double epsilon = 0.001);

	/**
	 * Local path.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @param epsilon The weight of paths of length 3.
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopLCP(int k, double epsilon = 0.001);

	/**
	 * Leicht-Holme-Newman index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllLHN();

	/**
	 * Leicht-Holme-Newman index.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 */
	void predLHN(std::vector<EdgeScore> &edgeScores);

	/**
	 * Leicht-Holme-Newman index.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopLHN(int k);

	/**
	 * Preferential attachment index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllPAT();

	/**
	 * Preferential attachment index.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 */
	void predPAT(std::vector<EdgeScore> &edgeScores);

	/**
	 * Preferential attachment index.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopPAT(int k);

	/**
	 * Resource allocation index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllRAL();

	/**
	 * Resource allocation index.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 */
	void predRAL(std::vector<EdgeScore> &edgeScores);

	/**
	 * Resource allocation index.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopRAL(int k);

	/**
	 * Random predictor.
	 * @param seed The random number generator seed.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllRND(long int seed = 0);

	/**
	 * Random predictor.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 * @param seed The random number generator seed.
	 */
	void predRND(std::vector<EdgeScore> &edgeScores, long int seed = 0);

	/**
	 * Random predictor.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @param seed The random number generator seed.
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopRND(int k, long int seed = 0);

	/**
	 * Salton index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllSAI();

	/**
	 * Salton index.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 */
	void predSAI(std::vector<EdgeScore> &edgeScores);

	/**
	 * Salton index.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopSAI(int k);

	/**
	 * Stochastic blocking model.
	 * @param maxIter Maximum number of iterations.
	 * @param seed Seed of the random number generator.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllSBM(int maxIter = 1000, long int seed = 0);

	/**
	 * Stochastic blocking model.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 * @param maxIter Maximum number of iterations.
	 * @param seed Seed of the random number generator.
	 */
	void predSBM(std::vector<EdgeScore> &edgeScores, int maxIter = 1000,
			long int seed = 0);

	/**
	 * Stochastic blocking model.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @param maxIter Maximum number of iterations.
	 * @param seed Seed of the random number generator.
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopSBM(int k, int maxIter = 1000, long int seed =
			0);

	/**
	 * Shortest path predictor.
	 * @param seed The random number generator seed.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllSHP(long int seed = 0);

	/**
	 * Shortest path predictor.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 * @param seed The random number generator seed.
	 */
	void predSHP(std::vector<EdgeScore> &edgeScores, long int seed = 0);

	/**
	 * Shortest path predictor.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @param seed The random number generator seed.
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopSHP(int k, long int seed = 0);

	/**
	 * Sorensen index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScore> predAllSOI();

	/**
	 * Sorensen index.
	 * @param edgeScores A input vector of negative edges. The score of each edge will be written in the member score.
	 */
	void predSOI(std::vector<EdgeScore> &edgeScores);

	/**
	 * Sorensen index.
	 * @param k Number of edges to be returned (the actual number may be smaller).
	 * @return The top k negative edge scores (the actual size may of the output be smaller than k).
	 */
	std::vector<EdgeScore> predTopSOI(int k);

	/**
	 * Adamic Adar predictor.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllADAByID();

	/**
	 * Common neighbors.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllCNEByID();

	/**
	 * Cannistraci resource allocation.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllCRAByID();

	/**
	 * Encoder-classifier link predictor.
	 * @param encoderName The name of the encoder. Possible values are: DPW (DeepWalk), HMSM (Hidden Metric Space Model), LVS (LargeVis), LEM (Laplacian Eigenmaps), LIN (LINE), LLE (Locally Linear Embedding), MFC (Matrix Factorization), and N2V (Node2Vec).
	 * @param classifierName The name of the classifier. Possible values are: FFN (feed-forward neural network withn default architecture), LSVM (linear SVM), LGR (logistic regression), NVB (naive Bayes). All classifiers except logistic regression requirte compilation with mlpack.
	 * @param dim The dimension of the embedding space. If set to zero, the default value is used (the default dimension depends on the ecnoder).
	 * @param posRatio Ratio of positive edges used in the training of the classifier.
	 * @param negRatio Ratio of negative edges used in the training of the classifier.
	 * @param seed Seed of the random number generator.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllECLByID(std::string encoderName = "N2V",
			std::string classifierName = "LGR", int dim = 0, double posRatio =
					1.0, double negRatio = 1.0, long int seed = 0);

	/**
	 * Encoder-similarity measure link predictor.
	 * @param encoderName The name of the encoder. Possible values are: DPW (DeepWalk), HMSM (Hidden Metric Space Model), LVS (LargeVis), LEM (Laplacian Eigenmaps), LIN (LINE), LLE (Locally Linear Embedding), MFC (Matrix Factorization), and N2V (Node2Vec).
	 * @param simMeasureName The name of the similarity measure. Possible values are: CSM (cosine similarity), DTP (dot product), L1 (L1 similarity), L2 (L2 similarity), PRS (Pearson similarity).
	 * @param dim The dimension of the embedding space. If set to zero, the default value is used (the default dimension depends on the ecnoder).
	 * @param seed Seed of the random number generator.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllESMByID(std::string encoderName = "N2V",
			std::string simMeasureName = "L2", int dim = 0, long int seed = 0);

	/**
	 * Fast blocking model.
	 * @param maxIter Maximum number of iterations.
	 * @param seed Seed of the random number generator.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllFBMByID(int maxIter = 50, long int seed = 0);

	/**
	 * Hub depromoted index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllHDIByID();

	/**
	 * Hub promoted index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllHPIByID();

	/**
	 * Hierarchical random graph.
	 * @param nbBeans Number of bins in edge statistics histogram.
	 * @param nbSamples Number of samples to take for predictions.
	 * @param seed Seed of the random number generator.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllHRGByID(int nbBeans = 25, int nbSamples =
			10000, long int seed = 0);

	/**
	 * Hypermap.
	 * @param m The parameter m (see the algorithm description).
	 * @param L The parameter L (see the algorithm description).
	 * @param gamma The power law exponent gamma (see the algorithm description).
	 * @param zeta The parameter zeta (see the algorithm description).
	 * @param T The parameter T (see the algorithm description).
	 * @param seed The random number generator seed.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllHYPByID(double m = 1.5, double L = 1,
			double gamma = 2.1, double zeta = 1, double T = 0.8, long int seed =
					0);

	/**
	 * Jackard index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllJIDByID();

	/**
	 * A scalable popularity similarity link predictor proposed in: "Kerrache, S., Alharbi,
	 * R. & Benhidour, H. A Scalable Similarity-Popularity Link Prediction Method. Sci Rep 10, 6394 (2020)".
	 * @param horizLim Horizon limit.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllKABByID(int horizLim = 2);

	/**
	 * Local path.
	 * @param epsilon The weight of paths of length 3.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllLCPByID(double epsilon = 0.001);

	/**
	 * Leicht-Holme-Newman index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllLHNByID();

	/**
	 * Preferential attachment index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllPATByID();

	/**
	 * Resource allocation index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllRALByID();

	/**
	 * Random predictor.
	 * @param seed The random number generator seed.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllRNDByID(long int seed = 0);

	/**
	 * Salton index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllSAIByID();

	/**
	 * Stochastic blocking model.
	 * @param maxIter Maximum number of iterations.
	 * @param seed Seed of the random number generator.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllSBMByID(int maxIter = 1000,
			long int seed = 0);

	/**
	 * Shortest path predictor.
	 * @param seed The random number generator seed.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllSHPByID(long int seed = 0);

	/**
	 * Sorensen index.
	 * @return The score of all negative edges.
	 */
	std::vector<EdgeScoreByID> predAllSOIByID();

#ifdef LINKPRED_WITH_OPENMP
	bool isParallel() const {
		return parallel;
	}

	void setParallel(bool parallel) {
		this->parallel = parallel;
	}
#endif
};

}
/* namespace Simp */
}
/* namespace LinkPred */

#endif /* EVALUATOR_HPP_ */
