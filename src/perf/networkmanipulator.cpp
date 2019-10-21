/*
 * This file is part of LinkPred.
 *
 * LinkPred: A high performance library for link prediction in complex networks.
 * Copyright (C) 2017  by Said Kerrache.
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

#include "linkpred/perf/networkmanipulator.hpp"
#include "linkpred/core/graphtraversal.hpp"
#include "linkpred/utils/utilities.hpp"
#include <iterator>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

namespace LinkPred {

#ifdef WITH_OPENMP
template<typename NetworkT> bool NetworkManipulator<NetworkT>::parallel = false;
#endif

template<typename NetworkT> void NetworkManipulator<NetworkT>::createTestData(
		NetworkCSP refNet, NetworkSP obsNet,
		std::shared_ptr<std::vector<EdgeType>> remLinks,
		std::shared_ptr<std::vector<EdgeType>> addLinks, double remRatio,
		double addRatio, bool keepConnected, long int seed) {

	std::size_t nbEdgesToRemove = std::rint(remRatio * refNet->getNbEdges());
	RandomGen rng(seed);

	// Add nodes
	for (auto it = refNet->nodesBegin(); it != refNet->nodesEnd(); ++it) {
		if (obsNet->addNode(it->second).first != it->first) {
			throw std::logic_error("IDs must be the same");
		}
	}

	if (keepConnected) {
		auto n = refNet->getNbNodes();
		auto m = refNet->getNbEdges();
		BFS<NetworkT> bfs(refNet);
		Collector<NetworkT> col;
		bfs.traverse(refNet->nodesBegin()->first, col);
		if (col.getVisited().size() != n) {
			throw std::domain_error("Network disconnected");
		}
		std::size_t nbEdgesToRemove = std::rint(remRatio * m);
		if (nbEdgesToRemove > m - n + 1) { // The size of a spanning tree is always nbNodes-1
			throw std::out_of_range(
					"The network will become disconnected after removing the requested number of edges.");
		}

		logger(logDebug1, "Getting a random spanning tree...")
		std::set<EdgeType> rstEdges;
		NetworkManipulator::rst(refNet, rng.getInt(),
				std::inserter(rstEdges, rstEdges.begin()));
		logger(logDebug1, "Done")

		// Select the rest of edges randomly
		std::vector<EdgeType> nonRstEdges;
		nonRstEdges.reserve(refNet->getNbEdges() - rstEdges.size());
		for (auto it = refNet->edgesBegin(); it != refNet->edgesEnd(); ++it) {
			if ((rstEdges.count(*it) == 0)
					&& (rstEdges.count(NetworkT::reverseEdge(*it)) == 0)) {
				nonRstEdges.push_back(*it);
			}
		}
		Utilities::selectRandomInPlace(nonRstEdges.begin(), nonRstEdges.end(),
				nbEdgesToRemove, rng.getInt());

		// Add spanning tree edges.
		for (auto it = rstEdges.begin(); it != rstEdges.end(); ++it) {
			obsNet->addEdge(NetworkT::start(*it), NetworkT::end(*it));
		}

		for (auto it = nonRstEdges.begin();
				it != nonRstEdges.begin() + nbEdgesToRemove; ++it) {
			remLinks->push_back(*it);
		}

		for (auto it = nonRstEdges.begin() + nbEdgesToRemove;
				it != nonRstEdges.end(); ++it) {
			obsNet->addEdge(NetworkT::start(*it), NetworkT::end(*it));
		}

	} else {
		// Select the rest of edges randomly
		std::vector<EdgeType> edgesCopy;
		edgesCopy.reserve(refNet->getNbEdges());
		for (auto it = refNet->edgesBegin(); it != refNet->edgesEnd(); ++it) {
			edgesCopy.push_back(*it);
		}
		Utilities::selectRandomInPlace(edgesCopy.begin(), edgesCopy.end(),
				nbEdgesToRemove, rng.getInt());

		for (auto it = edgesCopy.begin();
				it != edgesCopy.begin() + nbEdgesToRemove; ++it) {
			remLinks->push_back(*it);
		}
		for (auto it = edgesCopy.begin() + nbEdgesToRemove;
				it != edgesCopy.end(); ++it) {
			obsNet->addEdge(NetworkT::start(*it), NetworkT::end(*it));
		}
	}

	// Adding edges randomly
	for (auto it = refNet->rndNonEdgesBegin(addRatio, rng.getInt());
			it != refNet->rndNonEdgesEnd(); ++it) {
		obsNet->addEdge(NetworkT::start(*it), NetworkT::end(*it));
		addLinks->push_back(*it);
	}
	obsNet->assemble();
//	std::cout << "Edges:\n";
//	for (auto it = obsNet->edgesBegin(); it != obsNet->edgesEnd(); ++it) {
//		std::cout << obsNet->getLabel(obsNet->start(*it)) << "\t" << obsNet->getLabel(obsNet->end(*it)) << std::endl;
//	}
//	std::cout << "Removed links:\n";
//	for (auto it = remLinks->begin(); it != remLinks->end(); ++it) {
//		std::cout << obsNet->getLabel(obsNet->start(*it)) << "\t" << obsNet->getLabel(obsNet->end(*it)) << std::endl;
//	}
//	std::cout << "Nonedges:\n";
//	for (auto it = obsNet->nonEdgesBegin(); it != obsNet->nonEdgesEnd(); ++it) {
//		std::cout << obsNet->getLabel(obsNet->start(*it)) << "\t" << obsNet->getLabel(obsNet->end(*it)) << std::endl;
//	}
//	std::cout << "----------------------------------\n";
}

template<typename NetworkT> TestData<NetworkT,
		std::vector<typename NetworkManipulator<NetworkT>::EdgeType>> NetworkManipulator<
		NetworkT>::createTestDataSeqInter(NetworkCSP firstNet,
		NetworkCSP secondNet, bool aTP, double tpRatio, bool aTN,
		double tnRatio, LinkClass posClass, LinkClass negClass, long int seed) {

	auto refNet = std::make_shared<NetworkT>();
	auto obsNet = std::make_shared<NetworkT>();

	// Check for the common nodes
	for (auto it = secondNet->nodesBegin(); it != secondNet->nodesEnd(); ++it) {
		auto fit = firstNet->findLabel(it->second);
		if (fit != firstNet->labelsEnd()) {
			refNet->addNode(it->second);
			obsNet->addNode(it->second);
		}
	}

	auto remLinks = std::make_shared<std::vector<EdgeType>>();
	// Add edges to refNet
	for (auto it = secondNet->edgesBegin(); it != secondNet->edgesEnd(); ++it) {
		// Check if both ends exist in both networks
		auto srcLabel = secondNet->getLabel(secondNet->start(*it));
		auto endLabel = secondNet->getLabel(secondNet->end(*it));
		auto fitSrc = refNet->findLabel(srcLabel);
		if (fitSrc == refNet->labelsEnd()) {
			continue;
		}
		auto fitEnd = refNet->findLabel(endLabel);
		if (fitEnd == refNet->labelsEnd()) {
			continue;
		}
		refNet->addEdge(fitSrc->second, fitEnd->second);
		if (!firstNet->isEdge(firstNet->getID(srcLabel),
				firstNet->getID(endLabel))) {
			remLinks->push_back(
					refNet->makeEdge(fitSrc->second, fitEnd->second));
		}
	}
	refNet->assemble();

	auto addLinks = std::make_shared<std::vector<EdgeType>>();
	// Add edges to obsNet
	for (auto it = firstNet->edgesBegin(); it != firstNet->edgesEnd(); ++it) {
		// Check if both ends exist in both networks
		auto srcLabel = firstNet->getLabel(firstNet->start(*it));
		auto endLabel = firstNet->getLabel(firstNet->end(*it));
		auto fitSrc = obsNet->findLabel(srcLabel);
		if (fitSrc == obsNet->labelsEnd()) {
			continue;
		}
		auto fitEnd = obsNet->findLabel(endLabel);
		if (fitEnd == obsNet->labelsEnd()) {
			continue;
		}
		obsNet->addEdge(fitSrc->second, fitEnd->second);
		if (!secondNet->isEdge(secondNet->getID(srcLabel),
				secondNet->getID(endLabel))) {
			addLinks->push_back(
					obsNet->makeEdge(fitSrc->second, fitEnd->second));
		}
	}
	obsNet->assemble();

	RandomGen rng(seed);

	// The true positive links are the edges of obsNet that were not added
	std::set<EdgeType> addLinksSet;
	for (auto it = addLinks->begin(); it != addLinks->end(); ++it) {
		addLinksSet.insert(*it);
	}
	auto tpLinks = std::make_shared<std::vector<EdgeType>>();
	if (aTP) {
		for (auto it = obsNet->edgesBegin(); it != obsNet->edgesEnd(); ++it) {
			if ((addLinksSet.count(*it) == 0)
					&& (addLinksSet.count(NetworkT::reverseEdge(*it)) == 0)) {
				tpLinks->push_back(*it);
			}
		}
	} else {
		if ((tpRatio < 0) || (tpRatio > 1)) {
			logger(logError, "tpRatio: " << tpRatio)
			throw std::domain_error("The ratios must be between 0 and 1");
		}

		for (auto it = obsNet->rndEdgesBegin(tpRatio, rng.getInt());
				it != obsNet->rndEdgesEnd(); ++it) {
			if ((addLinksSet.count(*it) == 0)
					&& (addLinksSet.count(NetworkT::reverseEdge(*it)) == 0)) {
				tpLinks->push_back(*it);
			}
		}
	}

	// The true negative links are the non-edges of refNet that were not added
	auto tnLinks = std::make_shared<std::vector<EdgeType>>();
	if (aTN) {
		for (auto it = refNet->nonEdgesBegin(); it != refNet->nonEdgesEnd();
				++it) {
			if ((addLinksSet.count(*it) == 0)
					&& (addLinksSet.count(NetworkT::reverseEdge(*it)) == 0)) {
				tnLinks->push_back(*it);
			}
		}
	} else {
		if ((tnRatio < 0) || (tnRatio > 1)) {
			logger(logError, "tnRatio: " << tnRatio)
			throw std::domain_error("The ratios must be between 0 and 1");
		}

		for (auto it = refNet->rndNonEdgesBegin(tnRatio, rng.getInt());
				it != refNet->rndNonEdgesEnd(); ++it) {
			if ((addLinksSet.count(*it) == 0)
					&& (addLinksSet.count(NetworkT::reverseEdge(*it)) == 0)) {
				tnLinks->push_back(*it);
			}
		}
	}

	TestData<NetworkT, std::vector<EdgeType>> testData(refNet, obsNet, remLinks,
			addLinks, tpLinks, tnLinks, posClass, negClass);
	return testData;
}

//template<typename NetworkT> TestData<NetworkT,
//		std::vector<typename NetworkManipulator<NetworkT>::EdgeType>> NetworkManipulator<
//		NetworkT>::createTestDataSeq(NetworkCSP firstNet, NetworkCSP secondNet,
//		bool aTP, double tpRatio, bool aTN, double tnRatio, LinkClass posClass,
//		LinkClass negClass, long int seed) {
//
//	auto refNet = std::make_shared<NetworkT>();
//	auto obsNet = firstNet; // firstNet is the training network
//
//	// Now we create refNet
//	for (auto it = firstNet->nodesBegin(); it != firstNet->nodesEnd(); ++it) {
//		refNet->addNode(it->second); // Must have the same internal IDs
//	}
//
//	auto remLinks = std::make_shared<std::vector<EdgeType>>();
//	// Add edges to refNet
//	for (auto it = secondNet->edgesBegin(); it != secondNet->edgesEnd(); ++it) {
//		// Check if both ends exist in refNet (also first network)
//		auto srcLabel = secondNet->getLabel(secondNet->start(*it));
//		auto endLabel = secondNet->getLabel(secondNet->end(*it));
//		auto fitSrc = refNet->findLabel(srcLabel);
//		if (fitSrc == refNet->labelsEnd()) {
//			continue;
//		}
//		auto fitEnd = refNet->findLabel(endLabel);
//		if (fitEnd == refNet->labelsEnd()) {
//			continue;
//		}
//		refNet->addEdge(fitSrc->second, fitEnd->second);
//		if (!firstNet->isEdge(fitSrc->second, fitEnd->second)) {
//			remLinks->push_back(
//					refNet->makeEdge(fitSrc->second, fitEnd->second));
//		}
//	}
//	refNet->assemble();
//
//	// Get the set of added links
//	auto addLinks = std::make_shared<std::vector<EdgeType>>();
//	for (auto it = obsNet->edgesBegin(); it != obsNet->edgesEnd(); ++it) {
//		if (!refNet->isEdge(*it)) {
//			addLinks->push_back(*it);
//		}
//	}
//
//	RandomGen rng(seed);
//
//	// The true positive links are the edges of obsNet that were not added
//	std::set<EdgeType> addLinksSet;
//	for (auto it = addLinks->begin(); it != addLinks->end(); ++it) {
//		addLinksSet.insert(*it);
//	}
//	auto tpLinks = std::make_shared<std::vector<EdgeType>>();
//	if (aTP) {
//		for (auto it = obsNet->edgesBegin(); it != obsNet->edgesEnd(); ++it) {
//			if ((addLinksSet.count(*it) == 0)
//					&& (addLinksSet.count(NetworkT::reverseEdge(*it)) == 0)) {
//				tpLinks->push_back(*it);
//			}
//		}
//	} else {
//		if ((tpRatio < 0) || (tpRatio > 1)) {
//			logger(logError, "tpRatio: " << tpRatio)
//			throw std::domain_error("The ratios must be between 0 and 1");
//		}
//
//		for (auto it = obsNet->rndEdgesBegin(tpRatio, rng.getInt());
//				it != obsNet->rndEdgesEnd(); ++it) {
//			if ((addLinksSet.count(*it) == 0)
//					&& (addLinksSet.count(NetworkT::reverseEdge(*it)) == 0)) {
//				tpLinks->push_back(*it);
//			}
//		}
//	}
//
//	// The true negative links are the non-edges of refNet that were not added
//	auto tnLinks = std::make_shared<std::vector<EdgeType>>();
//	if (aTN) {
//		for (auto it = refNet->nonEdgesBegin(); it != refNet->nonEdgesEnd();
//				++it) {
//			if ((addLinksSet.count(*it) == 0)
//					&& (addLinksSet.count(NetworkT::reverseEdge(*it)) == 0)) {
//				tnLinks->push_back(*it);
//			}
//		}
//	} else {
//		if ((tnRatio < 0) || (tnRatio > 1)) {
//			logger(logError, "tnRatio: " << tnRatio)
//			throw std::domain_error("The ratios must be between 0 and 1");
//		}
//
//		for (auto it = refNet->rndNonEdgesBegin(tnRatio, rng.getInt());
//				it != refNet->rndNonEdgesEnd(); ++it) {
//			if ((addLinksSet.count(*it) == 0)
//					&& (addLinksSet.count(NetworkT::reverseEdge(*it)) == 0)) {
//				tnLinks->push_back(*it);
//			}
//		}
//	}
//
//	TestData<NetworkT, std::vector<EdgeType>> testData(refNet, obsNet, remLinks,
//			addLinks, tpLinks, tnLinks, posClass, negClass);
//	return testData;
//}

template<typename NetworkT> TestData<NetworkT,
		std::vector<typename NetworkManipulator<NetworkT>::EdgeType>> NetworkManipulator<
		NetworkT>::createTestDataSeq(NetworkCSP firstNet, NetworkCSP secondNet,
		bool aTP, double tpRatio, bool aTN, double tnRatio, LinkClass posClass,
		LinkClass negClass, long int seed, bool preGenerateTPN) {

	auto refNet = std::make_shared<NetworkT>();
	auto obsNet = firstNet; // firstNet is the training network

	// Now we create refNet
	for (auto it = firstNet->nodesBegin(); it != firstNet->nodesEnd(); ++it) {
		refNet->addNode(it->second); // Must have the same internal IDs
	}

	auto remLinks = std::make_shared<std::vector<EdgeType>>();
	// Add edges to refNet
	for (auto it = secondNet->edgesBegin(); it != secondNet->edgesEnd(); ++it) {
		// Check if both ends exist in refNet (also first network)
		auto srcLabel = secondNet->getLabel(secondNet->start(*it));
		auto endLabel = secondNet->getLabel(secondNet->end(*it));
		auto fitSrc = refNet->findLabel(srcLabel);
		if (fitSrc == refNet->labelsEnd()) {
			continue;
		}
		auto fitEnd = refNet->findLabel(endLabel);
		if (fitEnd == refNet->labelsEnd()) {
			continue;
		}
		refNet->addEdge(fitSrc->second, fitEnd->second);
		if (!firstNet->isEdge(fitSrc->second, fitEnd->second)) {
			remLinks->push_back(
					refNet->makeEdge(fitSrc->second, fitEnd->second));
//			std::cout << "rem: " << srcLabel << "\t" << endLabel << std::endl;
		}
	}
	refNet->assemble();

	// Get the set of added links
	auto addLinks = std::make_shared<std::vector<EdgeType>>();
	for (auto it = obsNet->edgesBegin(); it != obsNet->edgesEnd(); ++it) {
		if (!refNet->isEdge(*it)) {
			addLinks->push_back(*it);
		}
	}

	RandomGen rng(seed);

	if (preGenerateTPN) {
		auto eg = std::make_shared<
				TestEdgeGenerator<NetworkT, std::vector<EdgeType>> >(refNet,
				obsNet, remLinks, addLinks, aTP, tpRatio, aTN, tnRatio,
				rng.getInt());

		auto tpLinks = std::make_shared<std::vector<EdgeType>>();
		auto tnLinks = std::make_shared<std::vector<EdgeType>>();
		eg->template generateTN(std::back_inserter(*tnLinks));
		TestData<NetworkT, std::vector<EdgeType>> testData(refNet, obsNet,
				remLinks, addLinks, tpLinks, tnLinks, posClass, negClass);

		logger(logDebug, "Done")
		return testData;
	} else {
		auto eg = std::make_shared<
				TestEdgeGenerator<NetworkT, std::vector<EdgeType>> >(refNet,
				obsNet, remLinks, addLinks, aTP, tpRatio, aTN, tnRatio,
				rng.getInt());

		TestData<NetworkT, std::vector<EdgeType>> testData(refNet, obsNet,
				remLinks, addLinks, eg, posClass, negClass);
		logger(logDebug, "Done")
		return testData;
	}
}

template<typename NetworkT> TestData<NetworkT,
		std::vector<typename NetworkManipulator<NetworkT>::EdgeType>> NetworkManipulator<
		NetworkT>::createTestDataRem(NetworkCSP refNet, double remRatio,
		bool keepConnected, bool aTP, double tpRatio, bool aTN, double tnRatio,
		long int seed, bool preGenerateTPN) {
	logger(logDebug, "Creating test data...")

	if ((remRatio < 0) || (remRatio > 1)) {
		throw std::domain_error("The ratios must be between 0 and 1");
	}

	RandomGen rng(seed);
	// Create new network
	auto obsNet = std::make_shared<NetworkT>();
	auto remLinks = std::make_shared<std::vector<EdgeType>>();
	auto addLinks = std::make_shared<std::vector<EdgeType>>();

	createTestData(refNet, obsNet, remLinks, addLinks, remRatio, 0,
			keepConnected, rng.getInt());

	// The true positive links are the edges of obsNet
	auto tpLinks = std::make_shared<std::vector<EdgeType>>();
	if (aTP) {
		for (auto it = obsNet->edgesBegin(); it != obsNet->edgesEnd(); ++it) {
			tpLinks->push_back(*it);
		}
	} else {
		if ((tpRatio < 0) || (tpRatio > 1)) {
			logger(logError, "tpRatio: " << tpRatio)
			throw std::domain_error("The ratios must be between 0 and 1");
		}

		for (auto it = obsNet->rndEdgesBegin(tpRatio, rng.getInt());
				it != obsNet->rndEdgesEnd(); ++it) {
			tpLinks->push_back(*it);
		}
	}

	if (preGenerateTPN) {
		// The true negative links are the non-edges of refNet
		auto tnLinks = std::make_shared<std::vector<EdgeType>>();
		if (aTN) {
			for (auto it = refNet->nonEdgesBegin(); it != refNet->nonEdgesEnd();
					++it) {
				tnLinks->push_back(*it);
			}
		} else {
			if ((tnRatio < 0) || (tnRatio > 1)) {
				logger(logError, "tnRatio: " << tnRatio)
				throw std::domain_error("The ratios must be between 0 and 1");
			}

			for (auto it = refNet->rndNonEdgesBegin(tnRatio, rng.getInt());
					it != refNet->rndNonEdgesEnd(); ++it) {
				tnLinks->push_back(*it);
			}
		}

		TestData<NetworkT, std::vector<EdgeType>> testData(refNet, obsNet,
				remLinks, addLinks, tpLinks, tnLinks, FN, TN);
		logger(logDebug, "Done")

		return testData;
	} else {

		auto eg = std::make_shared<
				TestEdgeGenerator<NetworkT, std::vector<EdgeType>> >(refNet,
				obsNet, remLinks, addLinks, aTP, tpRatio, aTN, tnRatio,
				rng.getInt());

		TestData<NetworkT, std::vector<EdgeType>> testData(refNet, obsNet,
				remLinks, addLinks, eg, FN, TN);
		logger(logDebug, "Done")
		return testData;
	}
}

template<typename NetworkT> TestData<NetworkT,
		std::vector<typename NetworkManipulator<NetworkT>::EdgeType>> NetworkManipulator<
		NetworkT>::createTestDataAdd(NetworkCSP refNet, double addRatio,
		bool aTP, double tpRatio, bool aTN, double tnRatio, long int seed,
		bool preGenerateTPN) {
	logger(logDebug, "Creating test data...")

	if ((addRatio < 0) || (addRatio > 1)) {
		throw std::domain_error("The ratios must be between 0 and 1");
	}

	RandomGen rng(seed);
	// Create new network
	auto obsNet = std::make_shared<NetworkT>();
	auto remLinks = std::make_shared<std::vector<EdgeType>>();
	auto addLinks = std::make_shared<std::vector<EdgeType>>();

	createTestData(refNet, obsNet, remLinks, addLinks, 0, addRatio, false,
			rng.getInt());

	if (preGenerateTPN) {
		// The true positive links are the edges of refNet
		auto tpLinks = std::make_shared<std::vector<EdgeType>>();
		if (aTP) {
			for (auto it = refNet->edgesBegin(); it != refNet->edgesEnd();
					++it) {
				tpLinks->push_back(*it);
			}
		} else {
			if ((tpRatio < 0) || (tpRatio > 1)) {
				throw std::domain_error("The ratios must be between 0 and 1");
			}

			for (auto it = refNet->rndEdgesBegin(tpRatio, rng.getInt());
					it != refNet->rndEdgesEnd(); ++it) {
				tpLinks->push_back(*it);
			}
		}

		// The true negative links are the non-edges of obsNet
		auto tnLinks = std::make_shared<std::vector<EdgeType>>();
		if (aTN) {
			for (auto it = obsNet->nonEdgesBegin(); it != obsNet->nonEdgesEnd();
					++it) {
				tnLinks->push_back(*it);
			}
		} else {
			if ((tnRatio < 0) || (tnRatio > 1)) {
				throw std::domain_error("The ratios must be between 0 and 1");
			}

			for (auto it = obsNet->rndNonEdgesBegin(tnRatio, rng.getInt());
					it != obsNet->rndNonEdgesEnd(); ++it) {
				tnLinks->push_back(*it);
			}
		}

		TestData<NetworkT, std::vector<EdgeType>> testData(refNet, obsNet,
				remLinks, addLinks, tpLinks, tnLinks, TP, FP);

		logger(logDebug, "Done")
		return testData;
	} else {
		auto eg = std::make_shared<
				TestEdgeGenerator<NetworkT, std::vector<EdgeType>> >(refNet,
				obsNet, remLinks, addLinks, aTP, tpRatio, aTN, tnRatio,
				rng.getInt());

		TestData<NetworkT, std::vector<EdgeType>> testData(refNet, obsNet,
				remLinks, addLinks, eg, TP, FP);
		logger(logDebug, "Done")
		return testData;
	}
}

template<typename NetworkT> TestData<NetworkT,
		std::vector<typename NetworkManipulator<NetworkT>::EdgeType>> NetworkManipulator<
		NetworkT>::createTestData(NetworkCSP refNet, double remRatio,
		double addRatio, bool keepConnected, bool aTP, double tpRatio, bool aTN,
		double tnRatio, LinkClass posClass, LinkClass negClass, long int seed,
		bool preGenerateTPN) {
	logger(logDebug, "Creating test data...")

	if ((remRatio < 0) || (remRatio > 1) || (addRatio < 0) || (addRatio > 1)) {
		throw std::domain_error("The ratios must be between 0 and 1");
	}

	RandomGen rng(seed);
	// Create new network
	auto obsNet = std::make_shared<NetworkT>();
	auto remLinks = std::make_shared<std::vector<EdgeType>>();
	auto addLinks = std::make_shared<std::vector<EdgeType>>();

	createTestData(refNet, obsNet, remLinks, addLinks, remRatio, addRatio,
			keepConnected, rng.getInt());

	if (preGenerateTPN) {
		auto eg = std::make_shared<
				TestEdgeGenerator<NetworkT, std::vector<EdgeType>> >(refNet,
				obsNet, remLinks, addLinks, aTP, tpRatio, aTN, tnRatio,
				rng.getInt());

		auto tpLinks = std::make_shared<std::vector<EdgeType>>();
		auto tnLinks = std::make_shared<std::vector<EdgeType>>();
		eg->template generateTN(std::back_inserter(*tnLinks));
		TestData<NetworkT, std::vector<EdgeType>> testData(refNet, obsNet,
				remLinks, addLinks, tpLinks, tnLinks, posClass, negClass);

		logger(logDebug, "Done")
		return testData;
	} else {
		auto eg = std::make_shared<
				TestEdgeGenerator<NetworkT, std::vector<EdgeType>> >(refNet,
				obsNet, remLinks, addLinks, aTP, tpRatio, aTN, tnRatio,
				rng.getInt());

		TestData<NetworkT, std::vector<EdgeType>> testData(refNet, obsNet,
				remLinks, addLinks, eg, posClass, negClass);
		logger(logDebug, "Done")
		return testData;
	}
}

template<typename NetworkT> std::pair<
		typename NetworkManipulator<NetworkT>::NetworkCSP,
		std::shared_ptr<
				std::vector<typename NetworkManipulator<NetworkT>::EdgeType> > > NetworkManipulator<
		NetworkT>::rndConExtract(NetworkCSP net, double ratio, long int seed) {
	logger(logDebug, "Randomly extracting edges...")
	if ((ratio < 0) || (ratio > 1)) {
		throw std::domain_error("The ratio must be between 0 and 1");
	}

	auto n = net->getNbNodes();
	auto m = net->getNbEdges();
	std::size_t nbEdgesToExtract = std::rint(ratio * m);
	if (nbEdgesToExtract > m - n + 1) { // The size of a spanning tree is always nbNodes-1
		throw std::out_of_range(
				"The network will become disconnected after removing the requested number of edges.");
	}

	RandomGen rng(seed);
	logger(logDebug1, "Getting a random spanning tree...")
// Get a random spanning tree
	std::set<EdgeType> rstEdges;
	NetworkManipulator::rst(net, rng.getInt(),
			std::inserter(rstEdges, rstEdges.begin()));
	logger(logDebug1, "Done")

// Select the rest of edges randomly
	std::vector<EdgeType> nonRstEdges;
	nonRstEdges.reserve(net->getNbEdges() - rstEdges.size());
	for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it) {
		if ((rstEdges.count(*it) == 0)
				&& (rstEdges.count(NetworkT::reverseEdge(*it)) == 0)) {
			nonRstEdges.push_back(*it);
		}
	}
	Utilities::selectRandomInPlace(nonRstEdges.begin(), nonRstEdges.end(),
			nbEdgesToExtract, rng.getInt());

// Create new network
	auto resNet = std::make_shared<NetworkT>();
// Add nodes
	for (auto it = net->nodesBegin(); it != net->nodesEnd(); ++it) {
		if (resNet->addNode(it->second).first != it->first) {
			throw std::logic_error("IDs must be the same");
		}
	}

// Add spanning tree edges.
	for (auto it = rstEdges.begin(); it != rstEdges.end(); ++it) {
		resNet->addEdge(NetworkT::start(*it), NetworkT::end(*it));
	}

	auto extractedEdges = std::make_shared<std::vector<EdgeType>>();
	extractedEdges->reserve(nbEdgesToExtract);
	for (auto it = nonRstEdges.begin();
			it != nonRstEdges.begin() + nbEdgesToExtract; ++it) {
		extractedEdges->push_back(*it);
	}
	for (auto it = nonRstEdges.begin() + nbEdgesToExtract;
			it != nonRstEdges.end(); ++it) {
		resNet->addEdge(NetworkT::start(*it), NetworkT::end(*it));
	}

	resNet->assemble();
	logger(logDebug, "Done")
	return std::make_pair(resNet, extractedEdges);
}

template<typename NetworkT> std::pair<
		typename NetworkManipulator<NetworkT>::NetworkCSP,
		std::shared_ptr<
				std::vector<typename NetworkManipulator<NetworkT>::EdgeType> > > NetworkManipulator<
		NetworkT>::rndExtract(NetworkCSP net, double ratio, long int seed) {
	logger(logDebug, "Randomly extracting edges...")
	if ((ratio < 0) || (ratio > 1)) {
		throw std::domain_error("The ratio must be between 0 and 1");
	}

	std::size_t nbEdgesToExtract = std::rint(ratio * net->getNbEdges());
	RandomGen rng(seed);
// Select the rest of edges randomly
	std::vector<EdgeType> edgesCopy;
	edgesCopy.reserve(net->getNbEdges());
	for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it) {
		edgesCopy.push_back(*it);
	}
	Utilities::selectRandomInPlace(edgesCopy.begin(), edgesCopy.end(),
			nbEdgesToExtract, rng.getInt());

// Create new network
	auto resNet = std::make_shared<NetworkT>();
// Add nodes
	for (auto it = net->nodesBegin(); it != net->nodesEnd(); ++it) {
		if (resNet->addNode(it->second).first != it->first) {
			throw std::logic_error("IDs must be the same");
		}
	}

	auto extractedEdges = std::make_shared<std::vector<EdgeType>>();
	extractedEdges->reserve(nbEdgesToExtract);
	for (auto it = edgesCopy.begin();
			it != edgesCopy.begin() + nbEdgesToExtract; ++it) {
		extractedEdges->push_back(*it);
	}
	for (auto it = edgesCopy.begin() + nbEdgesToExtract; it != edgesCopy.end();
			++it) {
		resNet->addEdge(NetworkT::start(*it), NetworkT::end(*it));
	}

	resNet->assemble();
	logger(logDebug, "Done")
	return std::make_pair(resNet, extractedEdges);
}

#define NETWORKMANIPULATOR_CPP
#include "linkpred/instantiations.hpp"
#undef NETWORKMANIPULATOR_CPP

}
/* namespace LinkPred */
