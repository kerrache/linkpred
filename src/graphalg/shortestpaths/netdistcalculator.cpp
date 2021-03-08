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
#include "linkpred/graphalg/shortestpaths/netdistcalculator.hpp"
#ifdef LINKPRED_WITH_OPENMP
#include <omp.h>
#endif

namespace LinkPred {

template<typename Network, typename DistType, typename NbHopsType> std::pair<
		DistType, NbHopsType> ESPDistCalculator<Network, DistType, NbHopsType>::getDist(
		NodeID const & i, NodeID const & j) {
	switch (cacheLevel) {
	case NoCache: {
		return dijkstra.getDist(i, lengthMapId,
				NetDistCalculator<Network, DistType, NbHopsType>::discDist,
				NetDistCalculator<Network, DistType, NbHopsType>::discNbHops)->at(
				j);
	}

	case NodeCache: {
		std::pair<DistType, NbHopsType> res;
		bool missed = false;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPDistCalculatorNodeCacheUpdate)
#endif
		{
			if (!nodeDistCache) {
				missed = true;
			} else if (cachedNodeId == i) {
				res = nodeDistCache->at(j);
			} else if (cachedNodeId == j) {
				res = nodeDistCache->at(i);
			} else {
				missed = true;
			}
		}

		if (missed) {
			auto dsts =
					dijkstra.getDist(i, lengthMapId,
							NetDistCalculator<Network, DistType, NbHopsType>::discDist,
							NetDistCalculator<Network, DistType, NbHopsType>::discNbHops);
			res = dsts->at(j);

#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPDistCalculatorNodeCacheUpdate)
#endif
			{
				nodeDistCache = dsts;
				cachedNodeId = i;
			}
		}
		return res;
	}

	case NetworkCache: {
		std::pair<DistType, NbHopsType> res;
		bool missed = false;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPDistCalculatorNetworkCacheUpdate)
#endif
		{
			auto fit = netDistCache.find(i);
			if (fit != netDistCache.end()) {
				res = fit->second->at(j);
			} else {
				fit = netDistCache.find(j);
				if (fit != netDistCache.end()) {
					res = fit->second->at(i);
				} else {
					missed = true;
				}
			}
		}

		if (missed) {
			auto dist =
					dijkstra.getDist(i, lengthMapId,
							NetDistCalculator<Network, DistType, NbHopsType>::discDist,
							NetDistCalculator<Network, DistType, NbHopsType>::discNbHops);

#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPDistCalculatorNetworkCacheUpdate)
#endif
			{
				netDistCache[i] = dist;
			}

			res = dist->at(j);
		}
		return res;
	}

	default: {
		logger(logError, "Unknown cache level. Abort")
		exit(1);
	}
	}
	return std::make_pair(-1, -1);
}

template<typename Network, typename DistType, typename NbHopsType> std::pair<
		DistType, NbHopsType> ESPDistCalculator<Network, DistType, NbHopsType>::getIndDist(
		NodeID const & i, NodeID const & j) {

	return dijkstra.getIndDist(i, j, lengthMapId,
			NetDistCalculator<Network, DistType, NbHopsType>::discDist,
			NetDistCalculator<Network, DistType, NbHopsType>::discNbHops);
}

// TODO Add OMP critical
template<typename Network, typename DistType, typename NbHopsType> typename ESPDistCalculator<
		Network, DistType, NbHopsType>::NodeDistMapSP ESPDistCalculator<
		Network, DistType, NbHopsType>::getDist(NodeID const & i) {

	switch (cacheLevel) {
	case NoCache: {
		return dijkstra.getDist(i, lengthMapId,
				NetDistCalculator<Network, DistType, NbHopsType>::discDist,
				NetDistCalculator<Network, DistType, NbHopsType>::discNbHops);
	}

	case NodeCache: {
		if (!nodeDistCache) {
			nodeDistCache =
					dijkstra.getDist(i, lengthMapId,
							NetDistCalculator<Network, DistType, NbHopsType>::discDist,
							NetDistCalculator<Network, DistType, NbHopsType>::discNbHops);
			cachedNodeId = i;
		}

		if (cachedNodeId == i) {
			return nodeDistCache;
		} else {
			nodeDistCache =
					dijkstra.getDist(i, lengthMapId,
							NetDistCalculator<Network, DistType, NbHopsType>::discDist,
							NetDistCalculator<Network, DistType, NbHopsType>::discNbHops);
			cachedNodeId = i;
			return nodeDistCache;
		}
	}

	case NetworkCache: {
		auto fit = netDistCache.find(i);
		if (fit != netDistCache.end()) {
			return fit->second;
		} else {
			auto dist =
					dijkstra.getDist(i, lengthMapId,
							NetDistCalculator<Network, DistType, NbHopsType>::discDist,
							NetDistCalculator<Network, DistType, NbHopsType>::discNbHops);
			netDistCache[i] = dist;
			return dist;
		}
	}

	default: {
		throw std::runtime_error("Unknown cache level");
	}
	}
	return nullptr;
}

template<typename Network, typename DistType, typename NbHopsType> typename ESPLDistCalculator<
		Network, DistType, NbHopsType>::NodeDistMapSP ESPLDistCalculator<
		Network, DistType, NbHopsType>::getDist(NodeID const & i) {

	auto sdst = getSDistMap(i);
	auto dst = net->template createNodeMapSP<std::pair<DistType, NbHopsType>>();
	for (auto it = net->nodesBegin(); it != net->nodesEnd(); ++it) {
		(*dst)[it->first] = sdst->at(it->first);
	}
	return dst;
}

template<typename Network, typename DistType, typename NbHopsType> std::pair<
		DistType, NbHopsType> ESPLDistCalculator<Network, DistType, NbHopsType>::getDist(
		NodeID const & i, NodeID const & j) {

	// Copy to local variables to allow for swapping
	NodeID il = i;
	NodeID jl = j;

	switch (cacheLevel) {
	case NoCache: {
		auto dst = dijkstra.getDistL(il, lengthMapId, lim, discDist,
				discNbHops);
		return dst->at(jl);
	}

	case NodeCache: {
		bool missed = false;
		NodeSDistMapSP dst;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLDistCalculatorNodeCacheUpdate)
#endif
		{
			auto fit = netDistCache.find(il);
			if (fit != netDistCache.end()) {
				dst = fit->second;
			} else {
				fit = netDistCache.find(jl);
				if (fit != netDistCache.end()) {
					dst = fit->second;
					auto k = il;
					il = jl;
					jl = k;
				} else {
					missed = true;
				}
			}
		}
		if (missed) {
			dst = dijkstra.getDistL(il, lengthMapId, lim, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLDistCalculatorNodeCacheUpdate)
#endif
			{
				if (cachedNodeIds.size() >= maxNbNodesInCache) { // Need to remove from queue
					auto k = cachedNodeIds.front();
					cachedNodeIds.pop();
					netDistCache.erase(k);
				}
				cachedNodeIds.push(il);
				netDistCache[il] = dst;
			}
		}
		return dst->at(jl);
	}

	case NetworkCache: {
		bool missed = false;
		NodeSDistMapSP dst;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLDistCalculatorNetworkCacheUpdate)
#endif
		{
			auto fit = netDistCache.find(il);
			if (fit != netDistCache.end()) {
				dst = fit->second;
			} else {
				fit = netDistCache.find(jl);
				if (fit != netDistCache.end()) {
					dst = fit->second;
					auto k = il;
					il = jl;
					jl = k;
				} else {
					missed = true;
				}
			}
		}

		if (missed) {
			dst = dijkstra.getDistL(il, lengthMapId, lim, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLDistCalculatorNetworkCacheUpdate)
#endif
			{
				netDistCache[il] = dst;
			}
		}

		return dst->at(jl);
	}

	default: {
		throw std::runtime_error("Unknown cache level. Abort");
	}
	}

	return std::make_pair(discDist, discNbHops); // Should never reach this
}

template<typename Network, typename DistType, typename NbHopsType> typename ESPLDistCalculator<
		Network, DistType, NbHopsType>::NodeSDistMapSP ESPLDistCalculator<
		Network, DistType, NbHopsType>::getSDistMap(NodeID const & i) {

	switch (cacheLevel) {
	case NoCache: {
		return dijkstra.getDistL(i, lengthMapId, lim, discDist, discNbHops);
	}

	case NodeCache: {
		bool missed = false;
		NodeSDistMapSP dst;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLDistCalculatorNodeCacheUpdate)
#endif
		{
			auto fit = netDistCache.find(i);
			if (fit != netDistCache.end()) {
				dst = fit->second;
			} else {
				missed = true;
			}
		}

		if (missed) {
			dst = dijkstra.getDistL(i, lengthMapId, lim, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLDistCalculatorNodeCacheUpdate)
#endif
			{
				if (cachedNodeIds.size() >= maxNbNodesInCache) { // Need to remove from queue
					auto k = cachedNodeIds.front();
					cachedNodeIds.pop();
					netDistCache.erase(k);
				}
				cachedNodeIds.push(i);
				netDistCache[i] = dst;
			}
		}
		return dst;
	}

	case NetworkCache: {
		bool missed = false;
		NodeSDistMapSP dst;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLDistCalculatorNetworkCacheUpdate)
#endif
		{
			auto fit = netDistCache.find(i);
			if (fit != netDistCache.end()) {
				dst = fit->second;
			} else {
				missed = true;
			}
		}

		if (missed) {
			dst = dijkstra.getDistL(i, lengthMapId, lim, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLDistCalculatorNetworkCacheUpdate)
#endif
			{
				netDistCache[i] = dst;
			}
		}

		return dst;
	}

	default: {
		throw std::runtime_error("Unknown cache level. Abort");
	}
	}
}

template<typename Network, typename DistType, typename NbHopsType> typename ESPLDistCalculator<
		Network, DistType, NbHopsType>::NodeSDistMapSP ESPLDistCalculator<
		Network, DistType, NbHopsType>::getFinDistMapNoNeighb(
		NodeID const & srcNode) {

	auto i = srcNode;
	auto defVal = std::make_pair(discDist, discNbHops);
	auto dstnn =
			net->template createNodeSMapSP<std::pair<DistType, NbHopsType>>(
					defVal);

// Obtain the sparse distance map of node i
	auto disti = getSDistMap(i);
// Go through all elements of this map and filter self and neighbors
	for (auto it = disti->begin(); it != disti->end(); ++it) {
		auto j = it->first;
		if (i != j && !net->isEdge(i, j)) {
			(*dstnn)[j] = it->second;
		}
	}
	return dstnn;
}

template<typename Network, typename DistType, typename NbHopsType> std::pair<
		DistType, NbHopsType> ASPDistCalculator<Network, DistType, NbHopsType>::getDist(
		NodeID const & i, NodeID const & j) {
	auto dst = NetDistCalculator<Network, DistType, NbHopsType>::discDist;
	auto nbHops = NetDistCalculator<Network, DistType, NbHopsType>::discNbHops;

	for (auto it = distToLandmarks.cbegin(); it != distToLandmarks.cend();
			++it) {
		auto dsts = it->second->at(i).first;
		auto nbHopss = it->second->at(i).second;
		auto dste = it->second->at(j).first;
		auto nbHopse = it->second->at(j).second;
		if (dsts + dste < dst) {
			dst = dsts + dste;
			nbHops = nbHopss + nbHopse;
		}
	}
	return std::make_pair(dst, nbHops);
}

template<typename Network, typename DistType, typename NbHopsType> typename ASPDistCalculator<
		Network, DistType, NbHopsType>::NodeDistMapSP ASPDistCalculator<
		Network, DistType, NbHopsType>::getDist(NodeID const & i) {
	auto net = dijkstra.getNet();
	auto dist =
			net->template createNodeMapSP<std::pair<DistType, NbHopsType>>();
	for (NodeID j = 0; j < net->getNbNodes(); j++) {
		(*dist)[j] = getDist(i, j);
	}
	return dist;
}

template<typename Network, typename DsimType, typename NbHopsType> std::pair<
		DsimType, NbHopsType> ESPDsimCalculator<Network, DsimType, NbHopsType>::getDist(
		NodeID const & i, NodeID const & j) {
	switch (cacheLevel) {
	case NoCache: {
		return dijkstra.getDsim(i, lengthMapId,
				NetDistCalculator<Network, DsimType, NbHopsType>::discDist,
				NetDistCalculator<Network, DsimType, NbHopsType>::discNbHops)->at(
				j);
	}

	case NodeCache: {
		std::pair<DsimType, NbHopsType> res;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPDsimCalculatorNodeCacheUpdate1)
#endif
		{
			if (!nodeDsimCache) {
				nodeDsimCache =
						dijkstra.getDsim(i, lengthMapId,
								NetDistCalculator<Network, DsimType, NbHopsType>::discDist,
								NetDistCalculator<Network, DsimType, NbHopsType>::discNbHops);
				cachedNodeId = i;
			}
			if (cachedNodeId == i) {
				res = nodeDsimCache->at(j);
			} else if (cachedNodeId == j) {
				res = nodeDsimCache->at(i);
			} else {
				nodeDsimCache =
						dijkstra.getDsim(i, lengthMapId,
								NetDistCalculator<Network, DsimType, NbHopsType>::discDist,
								NetDistCalculator<Network, DsimType, NbHopsType>::discNbHops);
				cachedNodeId = i;
			}
			res = nodeDsimCache->at(j);
		}
		return res;
	}

	case NetworkCache: {
		std::pair<DsimType, NbHopsType> res;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPDsimCalculatorNetworkCacheUpdate1)
#endif
		{
			auto fit = netDsimCache.find(i);
			if (fit != netDsimCache.end()) {
				res = fit->second->at(j);
			} else {
				fit = netDsimCache.find(j);
				if (fit != netDsimCache.end()) {
					res = fit->second->at(i);
				} else {
					auto dsim =
							dijkstra.getDsim(i, lengthMapId,
									NetDistCalculator<Network, DsimType,
											NbHopsType>::discDist,
									NetDistCalculator<Network, DsimType,
											NbHopsType>::discNbHops);
					netDsimCache[i] = dsim;
					res = dsim->at(j);
				}
			}
		}
		return res;
	}
	default: {
		logger(logError, "Unknown cache level. Abort")
		exit(1);
	}
	}
	return std::make_pair(-1, -1);
}

template<typename Network, typename DistType, typename NbHopsType> std::pair<
		DistType, NbHopsType> ESPDsimCalculator<Network, DistType, NbHopsType>::getIndDist(
		NodeID const & i, NodeID const & j) {

	return dijkstra.getIndDsim(i, j, lengthMapId,
			NetDistCalculator<Network, DistType, NbHopsType>::discDist,
			NetDistCalculator<Network, DistType, NbHopsType>::discNbHops);
}

// TODO add omp critical
template<typename Network, typename DsimType, typename NbHopsType> typename ESPDsimCalculator<
		Network, DsimType, NbHopsType>::NodeDsimMapSP ESPDsimCalculator<
		Network, DsimType, NbHopsType>::getDist(NodeID const & i) {
#ifdef LINKPRED_WITH_OPENMP
	int nb = omp_get_thread_num();
	if (nb > 1) {
		throw std::runtime_error("Not safe in parallel");
	}
#endif
	switch (cacheLevel) {
	case NoCache: {
		return dijkstra.getDsim(i, lengthMapId,
				NetDistCalculator<Network, DsimType, NbHopsType>::discDist,
				NetDistCalculator<Network, DsimType, NbHopsType>::discNbHops);
	}

	case NodeCache: {

//#ifdef LINKPRED_WITH_OPENMP
//#pragma omp critical (ESPDsimCalculatorNodeCacheUpdate2)
//#endif
//		{
		if (!nodeDsimCache) {
			nodeDsimCache =
					dijkstra.getDsim(i, lengthMapId,
							NetDistCalculator<Network, DsimType, NbHopsType>::discDist,
							NetDistCalculator<Network, DsimType, NbHopsType>::discNbHops);
			cachedNodeId = i;
		}

		if (cachedNodeId == i) {
			return nodeDsimCache;
		} else {
			nodeDsimCache =
					dijkstra.getDsim(i, lengthMapId,
							NetDistCalculator<Network, DsimType, NbHopsType>::discDist,
							NetDistCalculator<Network, DsimType, NbHopsType>::discNbHops);
			cachedNodeId = i;
			return nodeDsimCache;
		}
	}
//	}
	case NetworkCache: {
//#ifdef LINKPRED_WITH_OPENMP
//#pragma omp critical (ESPDsimCalculatorNetworkCacheUpdate2)
//#endif
//		{
		auto fit = netDsimCache.find(i);
		if (fit != netDsimCache.end()) {
			return fit->second;
		} else {
			auto dist =
					dijkstra.getDsim(i, lengthMapId,
							NetDistCalculator<Network, DsimType, NbHopsType>::discDist,
							NetDistCalculator<Network, DsimType, NbHopsType>::discNbHops);
			netDsimCache[i] = dist;
			return dist;
		}
	}
//	}
	default: {
		throw std::runtime_error("Unknown cache level");
	}
	}
	return nullptr;
}

template<typename Network, typename DsimType, typename NbHopsType> std::pair<
		DsimType, NbHopsType> ASPDsimCalculator<Network, DsimType, NbHopsType>::getDist(
		NodeID const & i, NodeID const & j) {
	auto dst = NetDistCalculator<Network, DsimType, NbHopsType>::discDist;
	auto nbHops = NetDistCalculator<Network, DsimType, NbHopsType>::discNbHops;

	for (auto it = dsimToLandmarks.cbegin(); it != dsimToLandmarks.cend();
			++it) {
		auto dsts = it->second->at(i).first;
		auto nbHopss = it->second->at(i).second;
		auto dste = it->second->at(j).first;
		auto nbHopse = it->second->at(j).second;
		if (dsts + dste < dst) {
			dst = dsts + dste;
			nbHops = nbHopss + nbHopse;
		}
	}
	return std::make_pair(dst, nbHops);
}

template<typename Network, typename DsimType, typename NbHopsType> typename ASPDsimCalculator<
		Network, DsimType, NbHopsType>::NodeDsimMapSP ASPDsimCalculator<
		Network, DsimType, NbHopsType>::getDist(NodeID const & i) {
	auto net = dijkstra.getNet();
	auto dist =
			net->template createNodeMapSP<std::pair<DsimType, NbHopsType>>();
	for (NodeID j = 0; j < net->getNbNodes(); j++) {
		(*dist)[j] = getDist(i, j);
	}
	return dist;
}

template<typename Network, typename DistType, typename NbHopsType> std::tuple<
		DistType, DistType, NbHopsType> ESPSimlCalculator<Network, DistType,
		NbHopsType>::getSiml(NodeID const & i, NodeID const & j,
		NodeDistMapSP const simi, NodeDistMapSP const simj) const {

	double simij = 1;
	auto dhij = simi->at(j);
	for (auto it = net->neighbBegin(j); it != net->neighbEnd(j); ++it) {
		auto res = simi->at(net->end(*it));
		auto h = res.second;
		if (h != discNbHops) {
			double hp = h + 1;
			double ds = res.first + (*length)[*it];
			ds = ds / avgEdgeLength;
			simij *= (1 - std::pow(lambda, ds * hp));
		}
	}
	simij = 1 - simij;

	double simji = 1;
	auto dhji = simj->at(i);
	for (auto it = net->neighbBegin(i); it != net->neighbEnd(i); ++it) {
		auto res = simj->at(net->end(*it));
		auto h = res.second;
		if (h != discNbHops) {
			double hp = h + 1;
			double ds = res.first + (*length)[*it];
			ds = ds / avgEdgeLength;
			simji *= (1 - std::pow(lambda, ds * hp));
		}
	}
	simji = 1 - simji;

	return std::make_tuple(0.5 * (simij + simji),
			0.5 * (dhij.first + dhji.first), (dhij.second + dhji.second) / 2);
}

template<typename Network, typename DistType, typename NbHopsType> std::tuple<
		DistType, DistType, NbHopsType> ESPSimlCalculator<Network, DistType,
		NbHopsType>::getDirSiml(NodeID const & i, NodeID const & j,
		NodeDistMapSP const simi) const {
	throw std::runtime_error("not implemented");
//	std::cout << "max edge length = " << maxEdgeLength << std::endl;
	double diam = std::log((double) net->getNbNodes());

	double simij = 1;
	auto dhij = simi->at(j);
	for (auto it = net->neighbBegin(j); it != net->neighbEnd(j); ++it) {
		auto res = simi->at(net->end(*it));
		auto h = res.second;
		if (h != discNbHops) {
			double hp = h + 1;
			double ds = res.first + (*length)[*it];
			ds = ds / avgEdgeLength;
			simij *= (1 - std::pow(lambda, ds * hp / diam)); //std::exp(1.0 * (-lambda * hp - (1 - lambda) * ds));
		}
	}
	simij = 1 - simij;
	return std::make_tuple(simij, dhij.first, dhij.second);
}

template<typename Network, typename DistType, typename NbHopsType> std::tuple<
		DistType, DistType, NbHopsType> ESPSimlCalculator<Network, DistType,
		NbHopsType>::getSiml(NodeID const & i, NodeID const & j) {

	switch (cacheLevel) {
	case NoCache: {
		auto simi = dijkstra.getDist(i, lengthMapId, discDist, discNbHops);
		auto simj = dijkstra.getDist(j, lengthMapId, discDist, discNbHops);
		auto res = getSiml(i, j, simi, simj);
		return res;
	}

	case NodeCache: {
		NodeDistMapSP simi, simj;
		bool missedi = false;
		bool missedj = false;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPSimlCalculatorNodeCacheUpdate)
#endif
		{
			auto fit = netSimlCache.find(i);
			if (fit != netSimlCache.end()) {
				simi = fit->second;
			} else {
				missedi = true;
			}

			fit = netSimlCache.find(j);
			if (fit != netSimlCache.end()) {
				simj = fit->second;
			} else {
				missedj = true;
			}
		}

		if (missedi) {
			simi = dijkstra.getDist(i, lengthMapId, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPSimlCalculatorNodeCacheUpdate)
#endif
			{
				if (cachedNodeIds.size() >= maxNbNodesInCache) { // Need to remove from queue
					auto k = cachedNodeIds.front();
					cachedNodeIds.pop();
					netSimlCache.erase(k);
				}
				cachedNodeIds.push(i);
				netSimlCache[i] = simi;
			}
		}

		if (missedj) {
			simj = dijkstra.getDist(j, lengthMapId, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPSimlCalculatorNodeCacheUpdate)
#endif
			{
				if (cachedNodeIds.size() >= maxNbNodesInCache) { // Need to remove from queue
					auto k = cachedNodeIds.front();
					cachedNodeIds.pop();
					netSimlCache.erase(k);
				}
				cachedNodeIds.push(j);
				netSimlCache[j] = simj;
			}
		}
		return getSiml(i, j, simi, simj);
	}

	case NetworkCache: {
		NodeDistMapSP simi, simj;
		bool missedi = false;
		bool missedj = false;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPSimlCalculatorNetworkCacheUpdate)
#endif
		{
			auto fit = netSimlCache.find(i);
			if (fit != netSimlCache.end()) {
				simi = fit->second;
			} else {
				missedi = true;
			}

			fit = netSimlCache.find(j);
			if (fit != netSimlCache.end()) {
				simj = fit->second;
			} else {
				missedj = true;
			}
		}

		if (missedi) {
			simi = dijkstra.getDist(i, lengthMapId, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPSimlCalculatorNetworkCacheUpdate)
#endif
			{
				netSimlCache[i] = simi;
			}
		}

		if (missedj) {
			simj = dijkstra.getDist(j, lengthMapId, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPSimlCalculatorNetworkCacheUpdate)
#endif
			{
				netSimlCache[j] = simj;
			}
		}

		return getSiml(i, j, simi, simj);
	}

	default: {
		logger(logError, "Unknown cache level. Abort")
		exit(1);
	}
	}

	return std::make_tuple(selfSiml, discDist, discNbHops);
}

template<typename Network, typename DistType, typename NbHopsType> std::tuple<
		DistType, DistType, NbHopsType> ESPSimlCalculator<Network, DistType,
		NbHopsType>::getDirSiml(NodeID const & i, NodeID const & j) {

	switch (cacheLevel) {
	case NoCache: {
		auto simi = dijkstra.getDist(i, lengthMapId, discDist, discNbHops);
		auto res = getDirSiml(i, j, simi);
		return res;
	}

	case NodeCache: {
		NodeDistMapSP simi;
		bool missedi = false;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPSimlCalculatorNodeCacheUpdate)
#endif
		{
			auto fit = netSimlCache.find(i);
			if (fit != netSimlCache.end()) {
				simi = fit->second;
			} else {
				missedi = true;
			}
		}

		if (missedi) {
			simi = dijkstra.getDist(i, lengthMapId, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPSimlCalculatorNodeCacheUpdate)
#endif
			{
				if (cachedNodeIds.size() >= maxNbNodesInCache) { // Need to remove from queue
					auto k = cachedNodeIds.front();
					cachedNodeIds.pop();
					netSimlCache.erase(k);
				}
				cachedNodeIds.push(i);
				netSimlCache[i] = simi;
			}
		}

		return getDirSiml(i, j, simi);
	}

	case NetworkCache: {
		NodeDistMapSP simi;
		bool missedi = false;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPSimlCalculatorNetworkCacheUpdate)
#endif
		{
			auto fit = netSimlCache.find(i);
			if (fit != netSimlCache.end()) {
				simi = fit->second;
			} else {
				missedi = true;
			}
		}

		if (missedi) {
			simi = dijkstra.getDist(i, lengthMapId, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPSimlCalculatorNetworkCacheUpdate)
#endif
			{
				netSimlCache[i] = simi;
			}
		}

		return getDirSiml(i, j, simi);
	}

	default: {
		logger(logError, "Unknown cache level. Abort")
		exit(1);
	}
	}

	return std::make_tuple(selfSiml, discDist, discNbHops);
}

/*
 *
 */
template<typename Network, typename DistType, typename NbHopsType> std::tuple<
		DistType, DistType, NbHopsType> ESPLSimlCalculator<Network, DistType,
		NbHopsType>::getSiml(NodeID const & i, NodeID const & j,
		NodeSDistMapSP const simi, NodeSDistMapSP const simj) const {

	double maxDeg = net->getMaxDeg();

	double simij = 1;
	auto dhij = simi->at(j);
	if (dhij.second == discNbHops) {
		return std::make_tuple(0, discDist, discNbHops); // Disconnected
	} else {
		for (auto it = net->neighbBegin(j); it != net->neighbEnd(j);
				++it) {
			auto res = simi->at(net->end(*it));
			auto h = res.second;
			if (h != discNbHops) {
				double theta;
				if (useHops) {
					theta = h + 1;
				} else {
					theta = res.first + (*length)[*it];
				}
				double k = net->getDeg(net->end(*it));
//				simij *= std::pow(k / maxDeg, theta);
				simij *= std::pow(phi(k) / phi(maxDeg), theta);
			}
		}
	}
	simij = 1 - simij;

	double simji = 1;
	auto dhji = simj->at(i);
	if (dhji.second == discNbHops) {
		return std::make_tuple(0, discDist, discNbHops); // Disconnected
	} else {
		for (auto it = net->neighbBegin(i); it != net->neighbEnd(i);
				++it) {
			auto res = simj->at(net->end(*it));
			auto h = res.second;
			if (h != discNbHops) {
				double theta;
				if (useHops) {
					theta = h + 1;
				} else {
					theta = res.first + (*length)[*it];
				}
				double k = net->getDeg(net->end(*it));
//				simji *= std::pow(k / maxDeg, theta);
				simji *= std::pow(phi(k) / phi(maxDeg), theta);
			}
		}
	}
	simji = 1 - simji;

	return std::make_tuple(0.5 * (simij + simji),
			0.5 * (dhij.first + dhji.first), (dhij.second + dhji.second) / 2);
}

template<typename Network, typename DistType, typename NbHopsType> std::tuple<
		DistType, DistType, NbHopsType> ESPLSimlCalculator<Network, DistType,
		NbHopsType>::getDirSiml(NodeID const & i, NodeID const & j,
		NodeSDistMapSP const simi) const {

	double simij = 1;
	auto dhij = simi->at(j);
	if (dhij.second == discNbHops) {
		return std::make_tuple(0, discDist, discNbHops); // Disconnected
	} else {

		for (auto it = net->neighbBegin(j); it != net->neighbEnd(j);
				++it) {
			auto res = simi->at(net->end(*it));
			auto h = res.second;
			if (h != discNbHops) {
				double hp = h + 1;
				double ds = res.first + (*length)[*it];
				ds = ds / avgEdgeLength;
				simij *= (1 - std::pow(lambda, ds * hp));
			}
		}
	}
	simij = 1 - simij;
	return std::make_tuple(simij, dhij.first, dhij.second);
}

template<typename Network, typename DistType, typename NbHopsType> std::tuple<
		DistType, DistType, NbHopsType> ESPLSimlCalculator<Network, DistType,
		NbHopsType>::getSiml(NodeID const & i, NodeID const & j) {

	switch (cacheLevel) {
	case NoCache: {
		auto simi = dijkstra.getDistL(i, lengthMapId, lim, discDist,
				discNbHops);
		auto simj = dijkstra.getDistL(j, lengthMapId, lim, discDist,
				discNbHops);
		auto res = getSiml(i, j, simi, simj);
		return res;
	}

	case NodeCache: {
		NodeSDistMapSP simi, simj;
		bool missedi = false;
		bool missedj = false;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLSimlCalculatorNodeCacheUpdate)
#endif
		{
			auto fit = netSimlCache.find(i);
			if (fit != netSimlCache.end()) {
				simi = fit->second;
			} else {
				missedi = true;
			}

			fit = netSimlCache.find(j);
			if (fit != netSimlCache.end()) {
				simj = fit->second;
			} else {
				missedj = true;
			}
		}

		if (missedi) {
			simi = dijkstra.getDistL(i, lengthMapId, lim, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLSimlCalculatorNodeCacheUpdate)
#endif
			{
				if (cachedNodeIds.size() >= maxNbNodesInCache) { // Need to remove from queue
					auto k = cachedNodeIds.front();
					cachedNodeIds.pop();
					netSimlCache.erase(k);
				}
				cachedNodeIds.push(i);
				netSimlCache[i] = simi;
			}
		}

		if (missedj) {
			simj = dijkstra.getDistL(j, lengthMapId, lim, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLSimlCalculatorNodeCacheUpdate)
#endif
			{
				if (cachedNodeIds.size() >= maxNbNodesInCache) { // Need to remove from queue
					auto k = cachedNodeIds.front();
					cachedNodeIds.pop();
					netSimlCache.erase(k);
				}
				cachedNodeIds.push(j);
				netSimlCache[j] = simj;
			}
		}
		return getSiml(i, j, simi, simj);
	}

	case NetworkCache: {
		NodeSDistMapSP simi, simj;
		bool missedi = false;
		bool missedj = false;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLSimlCalculatorNetworkCacheUpdate)
#endif
		{
			auto fit = netSimlCache.find(i);
			if (fit != netSimlCache.end()) {
				simi = fit->second;
			} else {
				missedi = true;
			}

			fit = netSimlCache.find(j);
			if (fit != netSimlCache.end()) {
				simj = fit->second;
			} else {
				missedj = true;
			}
		}

		if (missedi) {
			simi = dijkstra.getDistL(i, lengthMapId, lim, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLSimlCalculatorNetworkCacheUpdate)
#endif
			{
				netSimlCache[i] = simi;
			}
		}

		if (missedj) {
			simj = dijkstra.getDistL(j, lengthMapId, lim, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLSimlCalculatorNetworkCacheUpdate)
#endif
			{
				netSimlCache[j] = simj;
			}
		}

		return getSiml(i, j, simi, simj);
	}

	default: {
		logger(logError, "Unknown cache level. Abort")
		exit(1);
	}
	}

	return std::make_tuple(selfSiml, discDist, discNbHops);
}

template<typename Network, typename DistType, typename NbHopsType> typename ESPLSimlCalculator<
		Network, DistType, NbHopsType>::NodeSDistMapSP ESPLSimlCalculator<
		Network, DistType, NbHopsType>::getSDistMap(NodeID const & i) {

	switch (cacheLevel) {
	case NoCache: {
		return dijkstra.getDistL(i, lengthMapId, lim, discDist, discNbHops);
	}

	case NodeCache: {
		NodeSDistMapSP disti;
		bool missedi = false;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLSimlCalculatorNodeCacheUpdate)
#endif
		{
			auto fit = netSimlCache.find(i);
			if (fit != netSimlCache.end()) {
				disti = fit->second;
			} else {
				missedi = true;
			}
		}

		if (missedi) {
			disti = dijkstra.getDistL(i, lengthMapId, lim, discDist,
					discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLSimlCalculatorNodeCacheUpdate)
#endif
			{
				if (cachedNodeIds.size() >= maxNbNodesInCache) { // Need to remove from queue
					auto k = cachedNodeIds.front();
					cachedNodeIds.pop();
					netSimlCache.erase(k);
				}
				cachedNodeIds.push(i);
				netSimlCache[i] = disti;
			}
		}

		return disti;
	}

	case NetworkCache: {
		NodeSDistMapSP disti;
		bool missedi = false;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLSimlCalculatorNetworkCacheUpdate)
#endif
		{
			auto fit = netSimlCache.find(i);
			if (fit != netSimlCache.end()) {
				disti = fit->second;
			} else {
				missedi = true;
			}
		}

		if (missedi) {
			disti = dijkstra.getDistL(i, lengthMapId, lim, discDist,
					discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLSimlCalculatorNetworkCacheUpdate)
#endif
			{
				netSimlCache[i] = disti;
			}
		}

		return disti;
	}

	default: {
		logger(logError, "Unknown cache level. Abort")
		exit(1);
	}
	}

	return nullptr;
}

template<typename Network, typename DistType, typename NbHopsType> typename ESPLSimlCalculator<
		Network, DistType, NbHopsType>::NodeSDistMapSP ESPLSimlCalculator<
		Network, DistType, NbHopsType>::getNnzSimlMap(
		NodeID const & srcNode) {

	auto i = srcNode;
	auto defVal = std::make_pair(0.0, discNbHops);
	auto simlMap = net->template createNodeSMapSP<
			std::pair<DistType, NbHopsType>>(defVal);

// Obtain the sparse distance map of node i
	auto disti = getSDistMap(i);
// Go through all elements of this map and compute similarity
	for (auto it = disti->begin(); it != disti->end(); ++it) {
		auto j = it->first;
		auto fit = simlMap->find(j);
		if (fit == simlMap->end()) {
			auto distj = getSDistMap(j);
			auto simij = getSiml(i, j, disti, distj);
			(*simlMap)[j] = std::make_pair(std::get<0>(simij),
					std::get<2>(simij));
		}
	}
	return simlMap;
}

template<typename Network, typename DistType, typename NbHopsType> typename ESPLSimlCalculator<
		Network, DistType, NbHopsType>::NodeSDistMapSP ESPLSimlCalculator<
		Network, DistType, NbHopsType>::getNnzSimlMapNoNeighb(
		NodeID const & srcNode) {

	auto i = srcNode;
	auto defVal = std::make_pair(0.0, discNbHops);
	auto simlMap = net->template createNodeSMapSP<
			std::pair<DistType, NbHopsType>>(defVal);

// Obtain the sparse distance map of node i
	auto disti = getSDistMap(i);
// Go through all elements of this map and compute similarity
	for (auto it = disti->begin(); it != disti->end(); ++it) {
		auto j = it->first;
		if (i != j && !net->isEdge(i, j)) {
			auto fit = simlMap->find(j);
			if (fit == simlMap->end()) {
				auto distj = getSDistMap(j);
				auto simij = getSiml(i, j, disti, distj);
				(*simlMap)[j] = std::make_pair(std::get<0>(simij),
						std::get<2>(simij));
			}
		}
	}
	return simlMap;
}

template<typename Network, typename DistType, typename NbHopsType> std::tuple<
		DistType, DistType, NbHopsType> ESPLSimlCalculator<Network, DistType,
		NbHopsType>::getDirSiml(NodeID const & i, NodeID const & j) {

	switch (cacheLevel) {
	case NoCache: {
		auto simi = dijkstra.getDistL(i, lengthMapId, lim, discDist,
				discNbHops);
		auto res = getDirSiml(i, j, simi);
		return res;
	}

	case NodeCache: {
		NodeSDistMapSP simi;
		bool missedi = false;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLSimlCalculatorNodeCacheUpdate)
#endif
		{
			auto fit = netSimlCache.find(i);
			if (fit != netSimlCache.end()) {
				simi = fit->second;
			} else {
				missedi = true;
			}
		}

		if (missedi) {
			simi = dijkstra.getDistL(i, lengthMapId, lim, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLSimlCalculatorNodeCacheUpdate)
#endif
			{
				if (cachedNodeIds.size() >= maxNbNodesInCache) { // Need to remove from queue
					auto k = cachedNodeIds.front();
					cachedNodeIds.pop();
					netSimlCache.erase(k);
				}
				cachedNodeIds.push(i);
				netSimlCache[i] = simi;
			}
		}

		return getDirSiml(i, j, simi);
	}

	case NetworkCache: {
		NodeSDistMapSP simi;
		bool missedi = false;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLSimlCalculatorNetworkCacheUpdate)
#endif
		{
			auto fit = netSimlCache.find(i);
			if (fit != netSimlCache.end()) {
				simi = fit->second;
			} else {
				missedi = true;
			}
		}

		if (missedi) {
			simi = dijkstra.getDistL(i, lengthMapId, lim, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (ESPLSimlCalculatorNetworkCacheUpdate)
#endif
			{
				netSimlCache[i] = simi;
			}
		}

		return getDirSiml(i, j, simi);
	}

	default: {
		logger(logError, "Unknown cache level. Abort")
		exit(1);
	}
	}

	return std::make_tuple(selfSiml, discDist, discNbHops);
}

/*
 **********************************************************************
 **********************************************************************
 **********************************************************************
 **********************************************************************
 **********************************************************************
 **********************************************************************
 */

template<typename Network, typename DistType, typename NbHopsType> std::tuple<
		DistType, DistType, NbHopsType> ESPIndSimlCalculator<Network, DistType,
		NbHopsType>::getDirIndSiml(
		std::vector<std::pair<DistType, NbHopsType>> const & dikj) const {

	double diam = std::log((double) net->getNbNodes());

	double simij = 1;
	for (auto it = dikj.cbegin() + 1; it != dikj.cend(); ++it) {
		double ds = it->first / avgEdgeLength;
		double hp = it->second;
		simij *= (1 - std::pow(lambda, ds * hp / diam)); //std::pow(lambda, ds * hp / diam2); //std::exp(1.0 * (-lambda * hp - (1 - lambda) * ds));
	}
	simij = 1 - simij;
	return std::make_tuple(simij, dikj[0].first, dikj[0].second);
}

template<typename Network, typename DistType, typename NbHopsType> std::tuple<
		DistType, DistType, NbHopsType> ESPIndSimlCalculator<Network, DistType,
		NbHopsType>::getIndSiml(NodeID const & i, NodeID const & j) {

	std::tuple<DistType, DistType, NbHopsType> sij, sji;

	auto fit = dirSimlCache.find(std::make_pair(i, j));
	if (fit == dirSimlCache.end()) {
		auto simi = dijkstra.getIndDistToAll(i, j, lengthMapId, discDist,
				discNbHops);
		std::vector<std::pair<DistType, NbHopsType>> dikj;
		auto dij = simi->at(j);
		dikj.push_back(dij);
		for (auto it = net->neighbBegin(j); it != net->neighbEnd(j);
				++it) {
			auto res = simi->at(net->end(*it));
			dikj.push_back(
					std::make_pair(res.first + length->at(*it),
							res.second + 1));
		}
		dirSimlCache[std::make_pair(i, j)] = dikj;
		sij = getDirIndSiml(dikj);
	} else {
		sij = getDirIndSiml(fit->second);
	}

	fit = dirSimlCache.find(std::make_pair(j, i));
	if (fit == dirSimlCache.end()) {
		auto simj = dijkstra.getIndDistToAll(j, i, lengthMapId, discDist,
				discNbHops);
		std::vector<std::pair<DistType, NbHopsType>> djki;
		auto dji = simj->at(i);
		djki.push_back(dji);
		for (auto it = net->neighbBegin(i); it != net->neighbEnd(i);
				++it) {
			auto res = simj->at(net->end(*it));
			djki.push_back(
					std::make_pair(res.first + length->at(*it),
							res.second + 1));
		}
		dirSimlCache[std::make_pair(j, i)] = djki;
		sji = getDirIndSiml(djki);
	} else {
		sij = getDirIndSiml(fit->second);
	}

	return std::make_tuple(0.5 * (std::get<0>(sij) + std::get<0>(sji)),
			0.5 * (std::get<1>(sij) + std::get<1>(sji)),
			(std::get<2>(sij) + std::get<2>(sji)) / 2);
}

template<typename Network, typename DistType, typename NbHopsType> std::tuple<
		DistType, DistType, NbHopsType> ESPIndSimlCalculator<Network, DistType,
		NbHopsType>::getDirIndSiml(NodeID const & i, NodeID const & j) {

	std::tuple<DistType, DistType, NbHopsType> sij;

	auto fit = dirSimlCache.find(std::make_pair(i, j));
	if (fit == dirSimlCache.end()) {
		auto simi = dijkstra.getIndDistToAll(i, j, lengthMapId, discDist,
				discNbHops);
		std::vector<std::pair<DistType, NbHopsType>> dikj;
		auto dij = simi->at(j);
		dikj.push_back(dij);
		for (auto it = net->neighbBegin(j); it != net->neighbEnd(j);
				++it) {
			auto res = simi->at(net->end(*it));
			dikj.push_back(
					std::make_pair(res.first + length->at(*it),
							res.second + 1));
		}
		dirSimlCache[std::make_pair(i, j)] = dikj;
		sij = getDirIndSiml(dikj);
	} else {
		sij = getDirIndSiml(fit->second);
	}

	return sij;
}

////////////////////////////////////////////////////////////////////////////
template<typename Network, typename DistType, typename NbHopsType> typename DESPLDistCalculator<
		Network, DistType, NbHopsType>::NodeDistMapSP DESPLDistCalculator<
		Network, DistType, NbHopsType>::getDist(NodeID const & i) {

	auto sdst = getSDistMap(i);
	auto dst = net->template createNodeMapSP<std::pair<DistType, NbHopsType>>();
	for (auto it = net->nodesBegin(); it != net->nodesEnd(); ++it) {
		(*dst)[it->first] = sdst->at(it->first);
	}
	return dst;
}

template<typename Network, typename DistType, typename NbHopsType> std::pair<
		DistType, NbHopsType> DESPLDistCalculator<Network, DistType, NbHopsType>::getDist(
		NodeID const & i, NodeID const & j) {

	switch (cacheLevel) {
	case NoCache: {
		auto dst = dijkstra.getDistL(i, lengthMapId, lim, discDist, discNbHops);
		return dst->at(j);
	}

	case NodeCache: {
		bool missed = false;
		NodeSDistMapSP dst;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (DESPLDistCalculatorNodeCacheUpdate)
#endif
		{
			auto fit = netDistCache.find(i);
			if (fit != netDistCache.end()) {
				dst = fit->second;
			} else {
				missed = true;
			}
		}
		if (missed) {
			dst = dijkstra.getDistL(i, lengthMapId, lim, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (DESPLDistCalculatorNodeCacheUpdate)
#endif
			{
				if (cachedNodeIds.size() >= maxNbNodesInCache) { // Need to remove from queue
					auto k = cachedNodeIds.front();
					cachedNodeIds.pop();
					netDistCache.erase(k);
				}
				cachedNodeIds.push(i);
				netDistCache[i] = dst;
			}
		}
		return dst->at(j);
	}

	case NetworkCache: {
		bool missed = false;
		NodeSDistMapSP dst;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (DESPLDistCalculatorNetworkCacheUpdate)
#endif
		{
			auto fit = netDistCache.find(i);
			if (fit != netDistCache.end()) {
				dst = fit->second;
			} else {
				missed = true;
			}
		}

		if (missed) {
			dst = dijkstra.getDistL(i, lengthMapId, lim, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (DESPLDistCalculatorNetworkCacheUpdate)
#endif
			{
				netDistCache[i] = dst;
			}
		}

		return dst->at(j);
	}

	default: {
		throw std::runtime_error("Unknown cache level. Abort");
	}
	}

	return std::make_pair(discDist, discNbHops); // Should never reach this
}

template<typename Network, typename DistType, typename NbHopsType> typename DESPLDistCalculator<
		Network, DistType, NbHopsType>::NodeSDistMapSP DESPLDistCalculator<
		Network, DistType, NbHopsType>::getSDistMap(NodeID const & i) {

	switch (cacheLevel) {
	case NoCache: {
		return dijkstra.getDistL(i, lengthMapId, lim, discDist, discNbHops);
	}

	case NodeCache: {
		bool missed = false;
		NodeSDistMapSP dst;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (DESPLDistCalculatorNodeCacheUpdate)
#endif
		{
			auto fit = netDistCache.find(i);
			if (fit != netDistCache.end()) {
				dst = fit->second;
			} else {
				missed = true;
			}
		}

		if (missed) {
			dst = dijkstra.getDistL(i, lengthMapId, lim, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (DESPLDistCalculatorNodeCacheUpdate)
#endif
			{
				if (cachedNodeIds.size() >= maxNbNodesInCache) { // Need to remove from queue
					auto k = cachedNodeIds.front();
					cachedNodeIds.pop();
					netDistCache.erase(k);
				}
				cachedNodeIds.push(i);
				netDistCache[i] = dst;
			}
		}
		return dst;
	}

	case NetworkCache: {
		bool missed = false;
		NodeSDistMapSP dst;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (DESPLDistCalculatorNetworkCacheUpdate)
#endif
		{
			auto fit = netDistCache.find(i);
			if (fit != netDistCache.end()) {
				dst = fit->second;
			} else {
				missed = true;
			}
		}

		if (missed) {
			dst = dijkstra.getDistL(i, lengthMapId, lim, discDist, discNbHops);
#ifdef LINKPRED_WITH_OPENMP
#pragma omp critical (DESPLDistCalculatorNetworkCacheUpdate)
#endif
			{
				netDistCache[i] = dst;
			}
		}

		return dst;
	}

	default: {
		throw std::runtime_error("Unknown cache level. Abort");
	}
	}
}

template<typename Network, typename DistType, typename NbHopsType> typename DESPLDistCalculator<
		Network, DistType, NbHopsType>::NodeSDistMapSP DESPLDistCalculator<
		Network, DistType, NbHopsType>::getFinDistMapNoNeighb(
		NodeID const & srcNode) {

	auto i = srcNode;
	auto defVal = std::make_pair(discDist, discNbHops);
	auto dstnn =
			net->template createNodeSMapSP<std::pair<DistType, NbHopsType>>(
					defVal);

// Obtain the sparse distance map of node i
	auto disti = getSDistMap(i);
// Go through all elements of this map and filter self and neighbors
	for (auto it = disti->begin(); it != disti->end(); ++it) {
		auto j = it->first;
		if (i != j && !net->isEdge(i, j)) {
			(*dstnn)[j] = it->second;
		}
	}
	return dstnn;
}

#define NETDISTCALCULATOR_CPP
#include "linkpred/instantiations.hpp"
#undef NETDISTCALCULATOR_CPP

} /* namespace LinkPred */
