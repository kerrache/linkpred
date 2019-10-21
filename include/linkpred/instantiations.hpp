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
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * \file
 * @brief Contains explicit instantiations of template classes.
 */

// Core *********************************************
#ifdef UNETWORK_CPP
template class UNetwork<>;
template class UNetwork<unsigned int>;
#endif

#ifdef DNETWORK_CPP
template class DNetwork<>;
template class DNetwork<unsigned int>;
#endif

#ifdef GRAPHTRAVERSAL_CPP
template class Counter<>;
template class Collector<>;
template class BFS<>;
template class DFS<>;
template class Counter<UNetwork<unsigned int>>;
template class Collector<UNetwork<unsigned int>>;
template class BFS<UNetwork<unsigned int>>;
template class DFS<UNetwork<unsigned int>>;
template class Counter<DNetwork<>>;
template class Collector<DNetwork<>>;
template class BFS<DNetwork<>>;
template class DFS<DNetwork<>>;
template class Counter<DNetwork<unsigned int>>;
template class Collector<DNetwork<unsigned int>>;
template class BFS<DNetwork<unsigned int>>;
template class DFS<DNetwork<unsigned int>>;
#endif

#ifdef NETDISTCALCULATOR_CPP
template class ESPDistCalculator<UNetwork<>, double, std::size_t>;
template class ESPDistCalculator<UNetwork<>, int, std::size_t>;
template class ESPDistCalculator<UNetwork<>, long int, std::size_t>;
template class ESPDistCalculator<UNetwork<>, long long int, std::size_t>;
template class ASPDistCalculator<UNetwork<>, double, std::size_t>;
template class ASPDistCalculator<UNetwork<>, int, std::size_t>;
template class ASPDistCalculator<UNetwork<>, long int, std::size_t>;
template class ASPDistCalculator<UNetwork<>, long long int, std::size_t>;
template class ESPDsimCalculator<UNetwork<>, double, std::size_t>;
template class ESPDsimCalculator<UNetwork<>, int, std::size_t>;
template class ESPDsimCalculator<UNetwork<>, long int, std::size_t>;
template class ESPDsimCalculator<UNetwork<>, long long int, std::size_t>;
template class ASPDsimCalculator<UNetwork<>, double, std::size_t>;
template class ASPDsimCalculator<UNetwork<>, int, std::size_t>;
template class ASPDsimCalculator<UNetwork<>, long int, std::size_t>;
template class ASPDsimCalculator<UNetwork<>, long long int, std::size_t>;
template class ESPSimlCalculator<UNetwork<>, double, std::size_t>;
template class ESPLSimlCalculator<UNetwork<>, double, std::size_t>;
template class ESPIndSimlCalculator<UNetwork<>, double, std::size_t>;
template class ESPLDistCalculator<UNetwork<>, double, std::size_t>;
template class ESPLDistCalculator<UNetwork<>, int, std::size_t>;
template class ESPLDistCalculator<UNetwork<>, long int, std::size_t>;
template class ESPLDistCalculator<UNetwork<>, long long int, std::size_t>;
template class DESPLDistCalculator<DNetwork<>, double, std::size_t>;
template class DESPLDistCalculator<DNetwork<>, int, std::size_t>;
template class DESPLDistCalculator<DNetwork<>, long int, std::size_t>;
template class DESPLDistCalculator<DNetwork<>, long long int, std::size_t>;
#endif

#ifdef DIJKSTRA_CPP
template class Dijkstra<UNetwork<>, double, std::size_t>;
template class Dijkstra<UNetwork<>, int, std::size_t>;
template class Dijkstra<UNetwork<>, long int, std::size_t>;
template class Dijkstra<UNetwork<>, long long int, std::size_t>;
template class Dijkstra<DNetwork<>, double, std::size_t>;
template class Dijkstra<DNetwork<>, int, std::size_t>;
template class Dijkstra<DNetwork<>, long int, std::size_t>;
template class Dijkstra<DNetwork<>, long long int, std::size_t>;
#endif

#ifdef PERFEVALUATOR_CPP
template class PerfEvaluator<>;
template class PerfEvaluator<TestData<DNetwork<>>, DLPredictor<>>;
template class PerfEvalExp<>;
//template class PerfEvalExp<DNetwork<>, TestData<DNetwork<>>, DLPredictor<>>;
#endif

// Predictors ***************************************
#ifdef UADAPREDICTOR_CPP
template class UADAPredictor<>;
template class UADAPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef UCNEPREDICTOR_CPP
template class UCNEPredictor<>;
template class UCNEPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef UCRAPREDICTOR_CPP
template class UCRAPredictor<>;
template class UCRAPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef UCSTPREDICTOR_CPP
template class UCSTPredictor<>;
template class UCSTPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef UENSPREDICTOR_CPP
template class UENSPredictor<>;
template class UENSPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef USNSPREDICTOR_CPP
template class USNSPredictor<>;
template class USNSPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef UKABPREDICTOR_CPP
template class UKABPredictor<>;
template class UKABPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef UFBMPREDICTOR_CPP
template class UFBMPredictor<>;
template class UFBMPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef UHDIPREDICTOR_CPP
template class UHDIPredictor<>;
template class UHDIPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef UHPIPREDICTOR_CPP
template class UHPIPredictor<>;
template class UHPIPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef UHRGPREDICTOR_CPP
template class UHRGPredictor<>;
template class UHRGPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef UHYPPREDICTOR_CPP
template class UHYPPredictor<>;
template class UHYPPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef UJIDPREDICTOR_CPP
template class UJIDPredictor<>;
template class UJIDPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef ULCPPREDICTOR_CPP
template class ULCPPredictor<>;
template class ULCPPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef ULHNPREDICTOR_CPP
template class ULHNPredictor<>;
#endif

#ifdef UPATPREDICTOR_CPP
template class UPATPredictor<>;
template class UPATPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef URALPREDICTOR_CPP
template class URALPredictor<>;
template class URALPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef URNDPREDICTOR_CPP
template class URNDPredictor<>;
template class URNDPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef USAIPREDICTOR_CPP
template class USAIPredictor<>;
template class USAIPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef USBMPREDICTOR_CPP
template class USBMPredictor<>;
template class USBMPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef USHPPREDICTOR_CPP
template class USHPPredictor<>;
template class USHPPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef USOIPREDICTOR_CPP
template class USOIPredictor<>;
template class USOIPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef USUMPREDICTOR_CPP
template class USUMPredictor<>;
template class USUMPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef UPOPPREDICTOR_CPP
template class UPOPPredictor<>;
template class UPOPPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef UNEDPREDICTOR_CPP
template class UNEDPredictor<>;
template class UNEDPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef UPNDPREDICTOR_CPP
template class UPNDPredictor<>;
template class UPNDPredictor<UNetwork<>, typename UNetwork<>::NonEdgeIterator>;
#endif

#ifdef DCNEPREDICTOR_CPP
template class DCNEPredictor<>;
template class DCNEPredictor<DNetwork<>, typename DNetwork<>::NonEdgeIterator>;
#endif

#ifdef DADAPREDICTOR_CPP
template class DADAPredictor<>;
template class DADAPredictor<DNetwork<>, typename DNetwork<>::NonEdgeIterator>;
#endif

#ifdef DSAIPREDICTOR_CPP
template class DSAIPredictor<>;
template class DSAIPredictor<DNetwork<>, typename DNetwork<>::NonEdgeIterator>;
#endif

#ifdef DJIDPREDICTOR_CPP
template class DJIDPredictor<>;
template class DJIDPredictor<DNetwork<>, typename DNetwork<>::NonEdgeIterator>;
#endif

#ifdef DSOIPREDICTOR_CPP
template class DSOIPredictor<>;
template class DSOIPredictor<DNetwork<>, typename DNetwork<>::NonEdgeIterator>;
#endif

#ifdef DHDIPREDICTOR_CPP
template class DHDIPredictor<>;
template class DHDIPredictor<DNetwork<>, typename DNetwork<>::NonEdgeIterator>;
#endif

#ifdef DHPIPREDICTOR_CPP
template class DHPIPredictor<>;
template class DHPIPredictor<DNetwork<>, typename DNetwork<>::NonEdgeIterator>;
#endif

#ifdef DPATPREDICTOR_CPP
template class DPATPredictor<>;
template class DPATPredictor<DNetwork<>, typename DNetwork<>::NonEdgeIterator>;
#endif

#ifdef DLCPPREDICTOR_CPP
template class DLCPPredictor<>;
template class DLCPPredictor<DNetwork<>, typename DNetwork<>::NonEdgeIterator>;
#endif

#ifdef DLHNPREDICTOR_CPP
template class DLHNPredictor<>;
template class DLHNPredictor<DNetwork<>, typename DNetwork<>::NonEdgeIterator>;
#endif

// Performance **************************************
#ifdef NETWORKMANIPULATOR_CPP
template class TestData<>;
template class NetworkManipulator<>;
template class TestData<UNetwork<unsigned int>>;
template class NetworkManipulator<UNetwork<unsigned int>>;
template class TestData<DNetwork<>>;
template class NetworkManipulator<DNetwork<>>;
template class TestData<DNetwork<unsigned int>>;
template class NetworkManipulator<DNetwork<unsigned int>>;
#endif

