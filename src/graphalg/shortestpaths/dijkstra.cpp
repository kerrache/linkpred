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

#include "LinkPredConfig.hpp"
#include "linkpred/graphalg/shortestpaths/dijkstra.hpp"
#include "linkpred/core/ds/bheap.hpp"
#include <algorithm>
#include <map>
#ifdef LINKPRED_WITH_OPENMP
#include <omp.h>
#endif

namespace LinkPred {

template<typename Network, typename DistType, typename NbHopsType> typename Dijkstra<
		Network, DistType, NbHopsType>::LengthMapIdType Dijkstra<Network,
		DistType, NbHopsType>::registerLengthMap(EdgeLengthMapSP length) {

	auto lengthMapId = getRndLengthMapId();
	auto fit = lengthMaps.find(lengthMapId);
	while (fit != lengthMaps.end()) {
		lengthMapId = getRndLengthMapId();
		fit = lengthMaps.find(lengthMapId);
	}
	lengthMaps[lengthMapId] = length;
	return lengthMapId;
}

template<typename Network, typename DistType, typename NbHopsType> void Dijkstra<
		Network, DistType, NbHopsType>::unregisterLengthMap(
		LengthMapIdType const & lengthMapId) {
	lengthMaps.erase(lengthMapId);
}

template<typename Network, typename DistType, typename NbHopsType> typename Dijkstra<
		Network, DistType, NbHopsType>::NodeDistMapSP Dijkstra<Network,
		DistType, NbHopsType>::getDsim(NodeID const & srcId,
		LengthMapIdType lengthMapId, DistType discDist,
		NbHopsType discNbHops) const {
	logger(logDebug, "Running Dijkstra...")
	auto length = lengthMaps.at(lengthMapId);
	auto dist =
			net->template createNodeMapSP<std::pair<DistType, NbHopsType>>();
	BHeap<NodeID, DistType> heap;
	struct RunData {
		bool valid = false;
		NodeID prev = 0;
		std::size_t nbHops;
		DistType dist;
		DistType smlr = 0;
	};
	std::vector<RunData> runData;
	runData.resize(net->getNbNodes());

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (NodeID i = 0; i < net->getNbNodes(); ++i) {
		runData[i].nbHops = discNbHops;
		runData[i].dist = discDist;
	}
	heap.push(srcId, 0);
	RunData rd { true, srcId, 0, 0, 0 };
	runData[srcId] = rd;

	while (heap.size() > 0) {
		auto s = heap.top();
		heap.pop();

		for (auto it = net->neighbBegin(s.first);
				it != net->neighbEnd(s.first); ++it) {
			auto neighb = Network::end(*it);
			auto el = length->at(Network::makeEdge(s.first, neighb));
			auto dst = s.second + el;
			auto nbHops = runData.at(s.first).nbHops + 1;
			if (!runData[neighb].valid) {
				runData[neighb].valid = true;
				runData[neighb].dist = dst;
				runData[neighb].nbHops = nbHops;
				runData[neighb].prev = s.first;
				runData[neighb].smlr = std::exp(-(double) nbHops); //
				heap.push(neighb, dst);
			} else {
				runData[neighb].smlr += std::exp(-(double) nbHops); //
				if (dst < runData[neighb].dist) {
					runData[neighb].dist = dst;
					runData[neighb].nbHops = nbHops;
					heap.increase(neighb, dst);
				}
			}
		}
	}

	logger(logDebug1, "Copying distances...")
	std::transform(runData.begin(), runData.end(), dist->begin(),
			[](RunData const & rd) {return std::make_pair( static_cast<DistType>(rd.smlr), rd.nbHops);}); //rd.deg
	logger(logDebug1, "Done")
	logger(logDebug, "Done")
	return dist;
}

template<typename Network, typename DistType, typename NbHopsType> std::pair<
		DistType, NbHopsType> Dijkstra<Network, DistType, NbHopsType>::getIndDsim(
		NodeID const & srcId, NodeID const & dstId,
		LengthMapIdType lengthMapId, DistType discDist,
		NbHopsType discNbHops) const {
	logger(logDebug, "Running Dijkstra...")
	auto length = lengthMaps.at(lengthMapId);
	BHeap<NodeID, DistType> heap;
	struct RunData {
		bool valid = false;
		NodeID prev = 0;
		std::size_t nbHops;
		DistType dist;
		DistType smlr = 0;
	};
	std::vector<RunData> runData;
	runData.resize(net->getNbNodes());

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (NodeID i = 0; i < net->getNbNodes(); ++i) {
		runData[i].nbHops = discNbHops;
		runData[i].dist = discDist;
	}
	heap.push(srcId, 0);
	RunData rd { true, srcId, 0, 0, 0 };
	runData[srcId] = rd;

	// First, we deal with the neighbors of the source node,
	// because that's where the removed edge may appear.
	for (auto it = net->neighbBegin(srcId); it != net->neighbEnd(srcId);
			++it) {
		auto neighb = Network::end(*it);
		if (neighb != dstId) {
			auto dst = length->at(Network::makeEdge(srcId, neighb));
			if (!runData[neighb].valid) {
				runData[neighb].valid = true;
				runData[neighb].dist = dst;
				runData[neighb].nbHops = 1;
				runData[neighb].prev = srcId;
				runData[neighb].smlr = std::exp(-1.0); // Number of hops is 1
				heap.push(neighb, dst);
			}
		}
	}
	while (heap.size() > 0) {
		auto s = heap.top();
		heap.pop();

		for (auto it = net->neighbBegin(s.first);
				it != net->neighbEnd(s.first); ++it) {
			auto neighb = Network::end(*it);
			auto el = length->at(Network::makeEdge(s.first, neighb));
			auto dst = s.second + el;
			auto nbHops = runData.at(s.first).nbHops + 1;
			if (!runData[neighb].valid) {
				runData[neighb].valid = true;
				runData[neighb].dist = dst;
				runData[neighb].nbHops = nbHops;
				runData[neighb].prev = s.first;
				runData[neighb].smlr = std::exp(-(double) nbHops); //
				heap.push(neighb, dst);
			} else {
				runData[neighb].smlr += std::exp(-(double) nbHops); //
				if (dst < runData[neighb].dist) {
					runData[neighb].dist = dst;
					runData[neighb].nbHops = nbHops;
					heap.increase(neighb, dst);
				}
			}
		}
	}

	logger(logDebug, "Done")
	return std::make_pair(static_cast<DistType>(runData[dstId].smlr),
			runData[dstId].nbHops);
}

template<typename Network, typename DistType, typename NbHopsType> typename Dijkstra<
		Network, DistType, NbHopsType>::NodeDistMapSP Dijkstra<Network,
		DistType, NbHopsType>::getDist(NodeID const & srcId,
		LengthMapIdType lengthMapId, DistType discDist,
		NbHopsType discNbHops) const {

	logger(logDebug, "Running Dijkstra...")

	auto length = lengthMaps.at(lengthMapId);
	BHeap<NodeID, DistType> heap;
	struct RunData {
		bool valid = false;
		NodeID prev = 0;
		std::size_t nbHops;
		DistType dist;
	};
	std::vector<RunData> runData;
	runData.resize(net->getNbNodes());
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = runData.begin(); it < runData.end(); ++it) {
		it->nbHops = discNbHops;
		it->dist = discDist;
	}
	heap.push(srcId, 0);
	RunData rd { true, srcId, 0, 0 };
	runData[srcId] = rd;

	while (heap.size() > 0) {
		auto s = heap.top();
		heap.pop();
		logger(logDebug1, "Popping: " << s.first << "\t" << s.second)
		for (auto it = net->neighbBegin(s.first);
				it != net->neighbEnd(s.first); ++it) {
			auto neighb = Network::end(*it);
			auto el = length->at(Network::makeEdge(s.first, neighb));
			auto dst = s.second + el;
			if (!runData[neighb].valid) {
				runData[neighb].valid = true;
				runData[neighb].dist = dst;
				runData[neighb].nbHops = runData.at(s.first).nbHops + 1;
				runData[neighb].prev = s.first;
				heap.push(neighb, dst);
			} else {
				if (dst < runData[neighb].dist) {
					runData[neighb].dist = dst;
					runData[neighb].nbHops = runData.at(s.first).nbHops + 1;
					heap.increase(neighb, dst);
				}
			}
		}
	}
	logger(logDebug1, "Copying distances...")
	auto dist =
			net->template createNodeMapSP<std::pair<DistType, NbHopsType>>();
	std::transform(runData.begin(), runData.end(), dist->begin(),
			[](RunData const & rd) {return std::make_pair(rd.dist,rd.nbHops);});
	logger(logDebug1, "Done")
	logger(logDebug, "Done")

	return dist;
}

template<typename Network, typename DistType, typename NbHopsType> typename Dijkstra<
		Network, DistType, NbHopsType>::NodeSDistMapSP Dijkstra<Network,
		DistType, NbHopsType>::getDistL(NodeID const & srcId,
		LengthMapIdType lengthMapId, std::size_t lim, DistType discDist,
		NbHopsType discNbHops) const {

	logger(logDebug, "Running Dijkstra...")

	auto length = lengthMaps.at(lengthMapId);
	BHeap<NodeID, DistType> heap;
	struct RunData {
		NodeID prev = 0;
		std::size_t nbHops;
		DistType dist;
	};
	std::map<NodeID, RunData> runData;

	heap.push(srcId, 0);
	RunData rd { srcId, 0, 0 };
	runData[srcId] = rd;

	while (heap.size() > 0) {
		auto s = heap.top();
		auto srd = runData.at(s.first);
		heap.pop();
		logger(logDebug1, "Popping: " << s.first << "\t" << s.second)
		if (srd.nbHops < lim) {
			for (auto it = net->neighbBegin(s.first);
					it != net->neighbEnd(s.first); ++it) {
				auto neighb = Network::end(*it);
				auto el = length->at(Network::makeEdge(s.first, neighb));
				auto dst = s.second + el;

				auto fit = runData.find(neighb);
				if (fit == runData.end()) {
					RunData rd;
					rd.dist = dst;
					rd.nbHops = srd.nbHops + 1;
					rd.prev = s.first;
					runData[neighb] = rd;
					heap.push(neighb, dst);
				} else {
					if (dst < fit->second.dist) {
						fit->second.dist = dst;
						fit->second.nbHops = srd.nbHops + 1;
						heap.increase(neighb, dst);
					}
				}
			}
		}
	}
	logger(logDebug1, "Copying distances...")
	auto defVal = std::make_pair(discDist, discNbHops);

	auto dist = net->template createNodeSMapSP<std::pair<DistType, NbHopsType>>(
			defVal);

	for (auto it = runData.begin(); it != runData.end(); ++it) {
		(*dist)[it->first] = std::make_pair(it->second.dist, it->second.nbHops);
	}
	logger(logDebug1, "Done")
	logger(logDebug, "Done")

	return dist;
}

template<typename Network, typename DistType, typename NbHopsType> std::pair<
		typename Dijkstra<Network, DistType, NbHopsType>::PathTypeSP, double> Dijkstra<
		Network, DistType, NbHopsType>::getShortestPath(
		NodeID const & srcId, NodeID const & dstId,
		LengthMapIdType lengthMapId, DistType discDist,
		NbHopsType discNbHops) const {

	logger(logDebug, "Running Dijkstra...")
	auto length = lengthMaps.at(lengthMapId);
	auto dist =
			net->template createNodeMapSP<std::pair<DistType, NbHopsType>>();
	BHeap<NodeID, DistType> heap;
	struct RunData {
		bool valid = false;
		NodeID prev = 0;
		std::size_t nbHops;
		DistType dist;
	};
	std::vector<RunData> runData;
	runData.resize(net->getNbNodes());
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = runData.begin(); it < runData.end(); ++it) {
		it->nbHops = discNbHops;
		it->dist = discDist;
	}
	heap.push(srcId, 0);
	RunData rd { true, srcId, 0, 0 };
	runData[srcId] = rd;

	while (heap.size() > 0) {
		auto s = heap.top();
		heap.pop();
		if (s.first == dstId) {
			break;
		}
		logger(logDebug1, "Popping: " << s.first << "\t" << s.second)
		for (auto it = net->neighbBegin(s.first);
				it != net->neighbEnd(s.first); ++it) {
			auto neighb = Network::end(*it);
			auto el = length->at(Network::makeEdge(s.first, neighb));
			auto dst = s.second + el;
			if (!runData[neighb].valid) {
				runData[neighb].valid = true;
				runData[neighb].dist = dst;
				runData[neighb].nbHops = runData.at(s.first).nbHops + 1;
				runData[neighb].prev = s.first;
				heap.push(neighb, dst);
			} else {
				if (dst < runData[neighb].dist) {
					runData[neighb].dist = dst;
					runData[neighb].nbHops = runData.at(s.first).nbHops + 1;
					heap.increase(neighb, dst);
				}
			}
		}
	}
	logger(logDebug1, "Creating path...")
	auto path = std::make_shared<std::vector<NodeID>>();
	if (runData[dstId].valid) {
		auto i = dstId;
		while (i != srcId) {
			path->push_back(i);
			i = runData[i].prev;
		}
		path->push_back(srcId);
		std::reverse(path->begin(), path->end());
	}
	logger(logDebug1, "Done")
	logger(logDebug, "Done")
	return std::make_pair(path, runData[dstId].dist);
}

template<typename Network, typename DistType, typename NbHopsType> std::pair<
		DistType, NbHopsType> Dijkstra<Network, DistType, NbHopsType>::getIndDist(
		NodeID const & srcId, NodeID const & dstId,
		LengthMapIdType lengthMapId, DistType discDist,
		NbHopsType discNbHops) const {

	logger(logDebug, "Running Dijkstra...")

	if (srcId == dstId) { // First, we get rid of this annoying special case
		return std::make_pair(0, 0);
	}

	auto length = lengthMaps.at(lengthMapId);
	BHeap<NodeID, DistType> heap;
	struct RunData {
		bool valid = false;
		NodeID prev = 0;
		std::size_t nbHops;
		DistType dist;
	};
	std::vector<RunData> runData;
	runData.resize(net->getNbNodes());

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = runData.begin(); it < runData.end(); ++it) {
		it->nbHops = discNbHops;
		it->dist = discDist;
	}
	RunData rd { true, srcId, 0, 0 };
	runData[srcId] = rd;

// First, we deal with the neighbors of the source node,
// because that's where the removed edge may appear.
	for (auto it = net->neighbBegin(srcId); it != net->neighbEnd(srcId);
			++it) {
		auto neighb = Network::end(*it);
		if (neighb != dstId) {
			auto dst = length->at(Network::makeEdge(srcId, neighb));
			if (!runData[neighb].valid) {
				runData[neighb].valid = true;
				runData[neighb].dist = dst;
				runData[neighb].nbHops = 1;
				runData[neighb].prev = srcId;
				heap.push(neighb, dst);
			}
		}
	}

	while (heap.size() > 0) {
		auto s = heap.top();
		if (s.first == dstId) { // If we reach the goal, we stop
			break;
		}
		heap.pop();
		logger(logDebug1, "Popping: " << s.first << "\t" << s.second)
		for (auto it = net->neighbBegin(s.first);
				it != net->neighbEnd(s.first); ++it) {
			auto neighb = Network::end(*it);
			auto el = length->at(Network::makeEdge(s.first, neighb));
			auto dst = s.second + el;
			if (!runData[neighb].valid) {
				runData[neighb].valid = true;
				runData[neighb].dist = dst;
				runData[neighb].nbHops = runData.at(s.first).nbHops + 1;
				runData[neighb].prev = s.first;
				heap.push(neighb, dst);
			} else {
				if (dst < runData[neighb].dist) {
					runData[neighb].dist = dst;
					runData[neighb].nbHops = runData.at(s.first).nbHops + 1;
					heap.increase(neighb, dst);
				}
			}
		}
	}

	logger(logDebug, "Done")
	return std::make_pair(runData[dstId].dist, runData[dstId].nbHops);
}

template<typename Network, typename DistType, typename NbHopsType> typename Dijkstra<
		Network, DistType, NbHopsType>::NodeDistMapSP Dijkstra<Network,
		DistType, NbHopsType>::getIndDistToAll(NodeID const & srcId,
		NodeID const & dstId, LengthMapIdType lengthMapId,
		DistType discDist, NbHopsType discNbHops) const {

	logger(logDebug, "Running Dijkstra...")

	auto length = lengthMaps.at(lengthMapId);
	BHeap<NodeID, DistType> heap;
	struct RunData {
		bool valid = false;
		NodeID prev = 0;
		std::size_t nbHops;
		DistType dist;
	};
	std::vector<RunData> runData;
	runData.resize(net->getNbNodes());

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = runData.begin(); it < runData.end(); ++it) {
		it->nbHops = discNbHops;
		it->dist = discDist;
	}
	RunData rd { true, srcId, 0, 0 };
	runData[srcId] = rd;

	heap.push(srcId, 0);

	auto eij = net->makeEdge(srcId, dstId);
	auto eji = net->makeEdge(dstId, srcId);

	while (heap.size() > 0) {
		auto s = heap.top();
		heap.pop();
		logger(logDebug1, "Popping: " << s.first << "\t" << s.second)
		for (auto it = net->neighbBegin(s.first);
				it != net->neighbEnd(s.first); ++it) {
			if (*it != eij && *it != eji) {
				auto neighb = Network::end(*it);
				auto el = length->at(*it);
				auto dst = s.second + el;
				if (!runData[neighb].valid) {
					runData[neighb].valid = true;
					runData[neighb].dist = dst;
					runData[neighb].nbHops = runData.at(s.first).nbHops + 1;
					runData[neighb].prev = s.first;
					heap.push(neighb, dst);
				} else {
					if (dst < runData[neighb].dist) {
						runData[neighb].dist = dst;
						runData[neighb].nbHops = runData.at(s.first).nbHops + 1;
						heap.increase(neighb, dst);
					}
				}
			}
		}
	}

	logger(logDebug1, "Copying distances...")
	auto dist =
			net->template createNodeMapSP<std::pair<DistType, NbHopsType>>();
	std::transform(runData.begin(), runData.end(), dist->begin(),
			[](RunData const & rd) {return std::make_pair(rd.dist,rd.nbHops);});
	logger(logDebug1, "Done")
	logger(logDebug, "Done")
	return dist;
}

#define DIJKSTRA_CPP
#include "linkpred/instantiations.hpp"
#undef DIJKSTRA_CPP

} /* namespace LinkPred */
