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
 * @brief Contains the definition of a class that simplifies the evaluation of link prediction algorithms.
 */

#ifndef SIMPEVALUATOR_HPP_
#define SIMPEVALUATOR_HPP_

#include "LinkPredConfig.hpp"
#include "linkpred/simp/perfres.hpp"
#include "linkpred/predictors/predictors.hpp"
#include "linkpred/perf/perf.hpp"
#include <memory>
#include <vector>
#include <map>
#include <set>

namespace LinkPred {
namespace Simp {

/**
 * @brief A class that simplifies the evaluation of link prediction algorithms.
 */
class Evaluator {
private:
	/**
	 * @brief Factory class.
	 */
	class Factory: public PEFactory<> {
	private:
		/**
		 * Create an encoder based on its name.
		 * @param name The name of the encoder.
		 * @param net The network.
		 * @param dim The dimension of the embedding.
		 * @param seed Seed of the random number generator.
		 * @return The encoder.
		 */
		std::shared_ptr<Encoder<>> createEncoder(std::string name,
				std::shared_ptr<const LinkPred::UNetwork<>> net, int dim,
				long int seed);

		/**
		 * Create a classifier based on its name.
		 * @param name The name of the classifier.
		 * @param dim The dimension of the input (edge embedding).
		 * @param seed Seed of the random number generator.
		 * @return The classifier.
		 */
		std::shared_ptr<Classifier<>> createClassifier(std::string name,
				int dim, long int seed);

		/**
		 * Create a similarity measure based on its name.
		 * @param name The name of the similarity measure.
		 * @param seed Seed of the random number generator.
		 * @return The similarity measure.
		 */
		std::shared_ptr<SimMeasure> createSimMeasure(std::string name);

	public:
#ifdef LINKPRED_WITH_OPENMP
		bool parallel = false; /**< Enable/disable parallelism. */
#endif
		std::vector<std::shared_ptr<ULPredictor<>>> predictors; /**< The predictors to be evaluated. */

		std::set<std::string> ada; /**< Evaluate ADA. */
		std::set<std::string> cne; /**< Evaluate CNE. */
		std::set<std::string> cra; /**< Evaluate CRA. */
		/**
		 * @brief Parameters of ECL.
		 */
		struct ECLParams {
			std::string encoderName = "N2V"; /**< The name of the encoder. Possible values are: DPW (DeepWalk), HMSM (Hidden Metric Space Model), LVS (LargeVis), LEM (Laplacian Eigenmaps), LIN (LINE), LLE (Locally Linear Embedding), MFC (Matrix Factorization) and N2V (Node2Vec). */
			std::string classifierName = "LGR"; /**< The name of the classifier. Possible values are: DTP (simple dot poroduct between the features), FFN (feed-forward neural network with default architecture), LSVM (linear SVM), LGR (logistic regression), NVB (naive Bayes). All classifiers except logistic regression requirte compilation with mlpack.*/
			int dim = 0; /**< The dimension of the embedding space. If set to zero, the default value is used (the default dimension depends on the encoder). */
			double posRatio = 1.0; /**< Ratio of positive edges used in the training of the classifier. */
			double negRatio = 1.0; /**< Ratio of negative edges used in the training of the classifier. */
			long int seed = 0; /**< Seed of the random number generator. */
		};
		std::map<std::string, ECLParams> ecl; /**< Evaluate ECL. */
		/**
		 * @brief Parameters of ESM.
		 */
		struct ESMParams {
			std::string encoderName = "N2V"; /**< The name of the encoder. Possible values are: DPW (DeepWalk), HMSM (Hidden Metric Space Model), LVS (LargeVis), LEM (Laplacian Eigenmaps), LIN (LINE), LLE (Locally Linear Embedding), MFC (Matrix Factorization) and N2V (Node2Vec). */
			std::string simMeasureName = "L2"; /**< The name of the similarity measure. Possible values are: CSM (cosine similarity), DTP (dot product), L1 (L1 similarity), L2 (L2 similarity), PRS (Pearson similarity). */
			int dim = 0; /**< The dimension of the embedding space. If set to zero, the default value is used (the default dimension depends on the encoder). */
			long int seed = 0; /**< Seed of the random number generator. */
		};
		std::map<std::string, ESMParams> esm; /**< Evaluate ESM. */
		/**
		 * @brief Parameters of FBM.
		 */
		struct FBMParams {
			int maxIter = 50; /**< Max iterations for FBM. */
			long int seed = 0; /**< Seed for FBM. */
		};
		std::map<std::string, FBMParams> fbm; /**< Evaluate FBM. */
		std::set<std::string> hdi; /**< Evaluate HDI. */
		std::set<std::string> hpi; /**< Evaluate HPI. */
		/**
		 * @brief Parameters of HRG.
		 */
		struct HRGParams {
			int nbBeans = 25; /**< Number of bins in edge statistics histogram in HRG. */
			int nbSamples = 10000; /**< Number of samples to take for predictions in HRG. */
			long int seed = 0; /**< Seed for HRG */
		};
		std::map<std::string, HRGParams> hrg; /**< Evaluate HRG. */
		/**
		 * @brief Parameters of HYP.
		 */
		struct HYPParams {
			double m = 1.5; /**< The parameter m (see the HYP algorithm description). */
			double L = 1; /**< The parameter L (see the HYP algorithm description). */
			double gamma = 2.1; /**< The power law exponent gamma (see the HYP algorithm description). */
			double zeta = 1; /**< The parameter zeta (see the HYP algorithm description). */
			double T = 0.8; /**< The parameter T (see the HYP algorithm description). */
			long int seed = 0; /**< The random number generator seed for HYP. */
		};
		std::map<std::string, HYPParams> hyp; /**< Evaluate HYP. */
		std::set<std::string> jid; /**< Evaluate JID. */
		/**
		 * @brief Parameters of KAB.
		 */
		struct KABParams {
			int horizLim = 2; /**< Horizon limit for KAB. */
		};
		std::map<std::string, KABParams> kab; /**< Evaluate KAB. */
		/**
		 * @brief Parameters of LCP.
		 */
		struct LCPParams {
			double epsilon = 0.001; /**< The weight of paths of length 3 in LCP. */
		};
		std::map<std::string, LCPParams> lcp; /**< Evaluate LCP. */
		std::set<std::string> lhn; /**< Evaluate LHN. */
		std::set<std::string> pat; /**< Evaluate PAT. */
		/**
		 * @brief Parameters of PST.
		 */
		struct PSTParams {
			std::string fileName; /**< File containing the scores of all non-exisityng links in the training network. */
		};
		std::map<std::string, PSTParams> pst; /**< Evaluate PST. */
		std::set<std::string> ral; /**< Evaluate RAL. */
		/**
		 * @brief Parameters of RND.
		 */
		struct RNDParams {
			long int seed = 0; /**< The random number generator seed for RND. */
		};
		std::map<std::string, RNDParams> rnd; /**< Evaluate RND. */
		std::set<std::string> sai; /**< Evaluate SAI. */
		/**
		 * @brief Parameters of SBM.
		 */
		struct SBMParams {
			int maxIter = 1000; /**< Max iterations for SBM. */
			long int seed = 0; /**< Seed for SBM. */
		};
		std::map<std::string, SBMParams> sbm; /**< Evaluate SBM. */
		/**
		 * @brief Parameters of SHP.
		 */
		struct SHPParams {
			long int seed = 0; /**< The random number generator seed for SHP. */
		};
		std::map<std::string, SHPParams> shp; /**< Evaluate SHP. */
		std::set<std::string> soi; /**< Evaluate SOI. */

		bool enableROC = false; /**< Enable ROC. */
		bool enablePR = false; /**< Enable PR. */
		bool enableTPR = false; /**< Enable TPR. */

		/**
		 * @param obsNet Observed network.
		 * @return The predictors.
		 */
		virtual std::vector<std::shared_ptr<ULPredictor<>>> getPredictors(
				std::shared_ptr<UNetwork<> const> obsNet);

		/**
		 * Creates performance measures.
		 * @param testData The test data.
		 * @return Perfomance measures.
		 */
		virtual std::vector<std::shared_ptr<PerfMeasure<>>> getPerfMeasures(
				TestData<> const &testData);

		/**
		 * Destructor.
		 */
		virtual ~Factory() = default;

	};

#ifdef LINKPRED_WITH_OPENMP
	bool parallel = false; /**< Enable/disable parallelism. */
#endif

	std::shared_ptr<Factory> factory; /**< The factory. */
	std::vector<std::vector<PerfRes>> perfRes; /**< The performance results. */

public:

	/**
	 * @param net The network.
	 */
	Evaluator();

	/**
	 * Copy constructor.
	 * @param that The object to copy.
	 */
	Evaluator(Evaluator const &that) = default;

	/**
	 * Destructor.
	 */
	virtual ~Evaluator() = default;

	/**
	 * Add Adamic Adar predictor.
	 * @param name A unique identifier for the predictor.
	 */
	void addADA(std::string const & name = "ADA");

	/**
	 * Add common neighbors.
	 * @param name A unique identifier for the predictor.
	 */
	void addCNE(std::string const & name = "CNE");

	/**
	 * Add Cannistraci resource allocation.
	 * @param name A unique identifier for the predictor.
	 */
	void addCRA(std::string const & name = "CRA");

	/**
	 * Add an encoder-classifier link predictor.
	 * @param name A unique identifier for the predictor.
	 * @param encoderName The name of the encoder. Possible values are: DPW (DeepWalk), HMSM (Hidden Metric Space Model), LVS (LargeVis), LEM (Laplacian Eigenmaps), LIN (LINE), LLE (Locally Linear Embedding), MFC (Matrix Factorization), and N2V (Node2Vec).
	 * @param classifierName The name of the classifier. Possible values are: FFN (feed-forward neural network withn default architecture), LSVM (linear SVM), LGR (logistic regression), NVB (naive Bayes). All classifiers except logistic regression requirte compilation with mlpack.
	 * @param dim The dimension of the embedding space. If set to zero, the default value is used (the default dimension depends on the ecnoder).
	 * @param posRatio Ratio of positive edges used in the training of the classifier.
	 * @param negRatio Ratio of negative edges used in the training of the classifier.
	 * @param seed Seed of the random number generator.
	 */
	void addECL(std::string const & name = "ECL-N2V-LGR",
			std::string encoderName = "N2V", std::string classifierName = "LGR",
			int dim = 0, double posRatio = 1.0, double negRatio = 1.0,
			long int seed = 0);

	/**
	 * Add an encoder-similarity measure link predictor.
	 * @param name A unique identifier for the predictor.
	 * @param encoderName The name of the encoder. Possible values are: DPW (DeepWalk), HMSM (Hidden Metric Space Model), LVS (LargeVis), LEM (Laplacian Eigenmaps), LIN (LINE), LLE (Locally Linear Embedding), MFC (Matrix Factorization), and N2V (Node2Vec).
	 * @param simMeasureName The name of the similarity measure. Possible values are: CSM (cosine similarity), DTP (dot product), L1 (L1 similarity), L2 (L2 similarity), PRS (Pearson similarity).
	 * @param dim The dimension of the embedding space. If set to zero, the default value is used (the default dimension depends on the ecnoder).
	 * @param seed Seed of the random number generator.
	 */
	void addESM(std::string const & name = "ESM-N2V-L2",
			std::string encoderName = "N2V", std::string simMeasureName = "L2",
			int dim = 0, long int seed = 0);

	/**
	 * Add fast blocking model.
	 * @param name A unique identifier for the predictor.
	 * @param maxIter Maximum number of iterations.
	 * @param seed Seed of the random number generator.
	 */
	void addFBM(std::string const & name = "FBM", int maxIter = 50,
			long int seed = 0);

	/**
	 * Add hub depromoted index.
	 * @param name A unique identifier for the predictor.
	 */
	void addHDI(std::string const & name = "HDI");

	/**
	 * Add Hub promoted index.
	 * @param name A unique identifier for the predictor.
	 */
	void addHPI(std::string const & name = "HPI");

	/**
	 * Add hierarchical random graph.
	 * @param name A unique identifier for the predictor.
	 * @param nbBeans Number of bins in edge statistics histogram.
	 * @param nbSamples Number of samples to take for predictions.
	 * @param seed Seed of the random number generator.
	 */
	void addHRG(std::string const & name = "HRG", int nbBeans = 25,
			int nbSamples = 10000, long int seed = 0);

	/**
	 * Add Hypermap.
	 * @param name A unique identifier for the predictor.
	 * @param m The parameter m (see the algorithm description).
	 * @param L The parameter L (see the algorithm description).
	 * @param gamma The power law exponent gamma (see the algorithm description).
	 * @param zeta The parameter zeta (see the algorithm description).
	 * @param T The parameter T (see the algorithm description).
	 * @param seed The random number generator seed.
	 */
	void addHYP(std::string const & name = "HYP", double m = 1.5, double L = 1,
			double gamma = 2.1, double zeta = 1, double T = 0.8, long int seed =
					0);

	/**
	 * Add Jackard index predictor.
	 * @param name A unique identifier for the predictor.
	 */
	void addJID(std::string const & name = "JID");

	/**
	 * Add the scalable popularity similarity link predictor proposed in: "Kerrache, S., Alharbi,
	 * R. & Benhidour, H. A Scalable Similarity-Popularity Link Prediction Method. Sci Rep 10, 6394 (2020)".
	 * @param name A unique identifier for the predictor.
	 * @param horizLim Horizon limit.
	 */
	void addKAB(std::string const & name = "KAB", int horizLim = 2);

	/**
	 * Add local path.
	 * @param name A unique identifier for the predictor.
	 * @param epsilon The weight of paths of length 3.
	 */
	void addLCP(std::string const & name = "LCP", double epsilon = 0.001);

	/**
	 * Add Leicht-Holme-Newman index.
	 * @param name A unique identifier for the predictor.
	 */
	void addLHN(std::string const & name = "LHN");

	/**
	 * Add preferential attachment index.
	 * @param name A unique identifier for the predictor.
	 */
	void addPAT(std::string const & name = "PAT");

	/**
	 * Add pre-stored results. This is designed to work with run(std::string const &obsEdgesFileName,
	 * std::string const &remEdgesFileName, long int seed) to measure the performance of user-defined
	 * prediction algorithms.
	 * @param name A unique identifier for the predictor.
	 * @param fileName File containing the scores of all non-exisityng links in the training network.
	 */
	void addPST(std::string const & name = "PST", std::string fileName =
			"pst.csv");

	/**
	 * Add resource allocation index.
	 * @param name A unique identifier for the predictor.
	 */
	void addRAL(std::string const & name = "RAL");

	/**
	 * Add random predictor.
	 * @param name A unique identifier for the predictor.
	 * @param seed The random number generator seed.
	 */
	void addRND(std::string const & name = "RND", long int seed = 0);

	/**
	 * Add Salton index.
	 * @param name A unique identifier for the predictor.
	 */
	void addSAI(std::string const & name = "SAI");

	/**
	 * Add stochastic blocking model.
	 * @param name A unique identifier for the predictor.
	 * @param maxIter Maximum number of iterations.
	 * @param seed Seed of the random number generator.
	 */
	void addSBM(std::string const & name = "SBM", int maxIter = 1000,
			long int seed = 0);

	/**
	 * Add shortest path predictor.
	 * @param name A unique identifier for the predictor.
	 * @param seed The random number generator seed.
	 */
	void addSHP(std::string const & name = "SHP", long int seed = 0);

	/**
	 * Add Sorensen index.
	 * @param name A unique identifier for the predictor.
	 */
	void addSOI(std::string const & name = "SOI");

	/**
	 * Add area under the ROC curve.
	 */
	void addROC();

	/**
	 * Add area under the precision-recall curve.
	 */
	void addPR();

	/**
	 * Add top precision.
	 */
	void addTPR();

	/**
	 * Generate test data and save it to file.
	 * @param fullNetFileName This the file containing the ground truth network.
	 * @param obsEdgesFileName The method writes the remaining edges into this file. This constitutes the training set.
	 * @param remEdgesFileName The method writes the removed edges into this file. This constitutes the positive examples of the test set.
	 * @param remRatio Edge remove ratio.
	 * @param keepConnected Whether to keep the graph connected when removing edges. This may be impossible for high ratios or if the network is initially disconnected.
	 * @param seed Seed for the andom number generator
	 */
	void genTestData(std::string const & fullNetFileName,
			std::string const &obsEdgesFileName,
			std::string const &remEdgesFileName, double remRatio = 0.1,
			bool keepConnected = false, long int seed = 0);

	/**
	 * Run a performance evaluation. This method creates test data by randomnly removing edges from the full network.
	 * @param fullNetFileName The network file name. This is the full network (ground truth) containing all edges.
	 * @param nbRuns Number of times the experiment is run. Each time the test set is changed.
	 * @param remRatio Edge remove ratio.
	 * @param keepConnected Whether to keep the graph connected when removing edges. This may be impossible for high ratios or if the network is initially disconnected.
	 * @param seed Seed for the andom number generator
	 */
	void run(std::string const &fullNetFileName, int nbRuns = 10,
			double remRatio = 0.1, bool keepConnected = false,
			long int seed = 0);

	/**
	 * Run a performance evaluation. This method uses the test data given as input
	 * @param obsEdgesFileName A file containing observed edges (training set).
	 * @param remEdgesFileName A file containing removed edges (the positive examples of the test set).
	 */
	void run(std::string const &obsEdgesFileName,
			std::string const &remEdgesFileName);

	/**
	 * @param iter Iteration number.
	 * @return The performance results at the specified iteration.
	 */
	std::vector<PerfRes> getPerfRes(int iter) {
		return perfRes[iter];
	}

#ifdef LINKPRED_WITH_OPENMP
	bool isParallel() const {
		return parallel;
	}

	void setParallel(bool parallel) {
		this->parallel = parallel;
		factory->parallel = parallel;
	}
#endif
};

}
/* namespace Simp */
}
/* namespace LinkPred */

#endif /* SIMPEVALUATOR_HPP_ */

