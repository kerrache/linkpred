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

/**
 * \file
 * @ingroup Utils
 * @brief Contains the implementation of miscellaneous useful methods.
 */

#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_

#include "LinkPredConfig.hpp"
#include "linkpred/utils/randomgen.hpp"
#include "linkpred/utils/log.hpp"
#include "linkpred/core/ds/lmapqueue.hpp"
#include <string>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <queue>
#include <fstream>
#include <iostream>
#include <string>
#include <limits>
#include <algorithm>
#include <iomanip>
#ifdef LINKPRED_WITH_MPI
#include <mpi.h>
#endif

namespace LinkPred {

/**
 * @brief Enumeration of different sorting orders.
 */
enum SortOrder {
	None, /**< Not sorted. */
	Inc, /**< Sorted in increasing order. */
	Dec /**< Sorted in decreasing order. */
};

/**
 * @brief Some utility functions.
 */
namespace Utils {

/**
 * Trim the string from left.
 * @param str Input string.
 * @param spaces List of spaces.
 */
inline std::string& ltrim(std::string& str,
		const char* spaces = " \t\n\r\f\v") {
	str.erase(0, str.find_first_not_of(spaces));
	return str;
}

/**
 * Trim the string from right.
 * @param str Input string.
 * @param spaces List of spaces.
 */
inline std::string& rtrim(std::string& str,
		const char* spaces = " \t\n\r\f\v") {
	str.erase(str.find_last_not_of(spaces) + 1);
	return str;
}

/**
 * Trim the string from left and right.
 * @param str Input string.
 * @param spaces List of spaces.
 */
inline std::string& trim(std::string& str, const char* spaces = " \t\n\r\f\v") {
	return ltrim(rtrim(str, spaces), spaces);
}

/**
 * Clear the queue (q becomes empty after the call).
 * @tparam T the data type stored in the queue.
 * @param q The queue.
 */
template<typename T> void clear(std::queue<T> & q) {
	std::queue<T> empty;
	std::swap(q, empty);
}

/**
 * Sort a range.
 * @tparam T The data type.
 * @tparam IteratorT The iterator type.
 * @param begin Iterator to the first element in the range.
 * @param end Iterator to one-past-the-last element in the range.
 * @param sortOrder The requested sorting order.
 */
template<typename IteratorT> void sort(IteratorT begin, IteratorT end,
		SortOrder sortOrder) {
	switch (sortOrder) {
	case Inc: {
		std::sort(begin, end);
		break;
	}
	case Dec: {
		std::sort(begin, end,
				std::greater<
						typename std::iterator_traits<IteratorT>::value_type>());
		break;
	}
	case None: {
		break;
	}
	default: {
		throw std::invalid_argument("Unknown sorting status");
	}
	}
}

/**
 * @param n Size of the permutation.
 * @return A vector containing a random permutation of [0..n-1].
 */
std::vector<std::size_t> getRndPerm(std::size_t n);

/**
 * @param n Size of the permutation.
 * @param seed The seed.
 * @return A vector containing a random permutation of [0..n-1].
 */
std::vector<std::size_t> getRndPerm(std::size_t n, long int seed);

/**
 * Fit a discrete power law.
 * @param data The data to be fitted.
 * @return A pair containing respectively, the power of the power law and the minimum value.
 */
std::pair<double, double> plFit(std::vector<std::size_t> const & data);

/**
 * Fit a continuous power law.
 * @param data The data to be fitted.
 * @return A pair containing respectively, the power of the power law and the minimum value.
 */
std::pair<double, double> plFit(std::vector<double> const & data);

/**
 * Print vector to standard output.
 * @param v The vector.
 * @param name String prefixed to the print.
 */
template<typename T> void print(std::vector<T> const & v, std::string name) {
	std::cout << name << ": " << std::endl;
	for (std::size_t i = 0; i < v.size(); i++) {
		std::cout << v[i] << "\t";
	}
	std::cout << std::endl;
}

/**
 * Flip a pair.
 * @param p The original pair.
 * @return The flipped pair.
 */
template<typename U, typename V> std::pair<V, U> flip(
		const std::pair<U, V> &p) {
	return std::pair<V, U>(p.second, p.first);
}

/**
 * Flip a map.
 * @param map The original map.
 * @return The flipped map.
 */
template<typename U, typename V> std::multimap<V, U> flipMap(
		const std::map<U, V> &map) {
	std::multimap<V, U> revMap;
	std::transform(map.begin(), map.end(),
			std::inserter(revMap, revMap.begin()), flip<U, V>);
	return revMap;
}

/**
 * Randomly select k elements without repetition using Floyd's selection algorithm.
 * @param begin Iterator pointing to the start of elements.
 * @param end Iterator pointing to the end of elements.
 * @param k The number of elements to select. If k is larger than the actual number of elements, std::out_of_range is thrown.
 * @param seed Seed for the random number generator.
 */
template<typename RandomIterator> std::set<
		typename std::iterator_traits<RandomIterator>::value_type> selectRandom(
		RandomIterator begin, RandomIterator end, std::size_t k,
		long int seed) {
	logger(logDebug, "Selecting elements randomly...")
	typedef typename std::iterator_traits<RandomIterator>::value_type T;
	std::size_t n = std::distance(begin, end);
	if (k > n) {
		throw std::out_of_range(
				"The number of elements to select is larger than the total number of elements.");
	}
	RandomGen rng(seed);
	std::set<T> selected;
	for (std::size_t i = n - k; i < n; i++) {
		std::size_t ind = rng.getSInt(0, i - 1);
		RandomIterator it = begin + ind;
		if (selected.count(*it) == 0) {
			selected.insert(*it);
		} else {
			selected.insert(*(begin + i));
		}
	}
	logger(logDebug, "Done")
	return selected;
}

/**
 * Randomly select k elements without repetition using Floyd's selection algorithm.
 * The selected elements are put at the end of the range.
 * @param begin Iterator pointing to the start of elements.
 * @param end Iterator pointing to the end of elements.
 * @param k The number of elements to select. If k is larger than the actual number of elements, std::out_of_range is thrown.
 * @param seed Seed for the random number generator.
 */
template<typename RandomIterator> void selectRandomInPlace(RandomIterator begin,
		RandomIterator end, std::size_t k, long int seed) {
	logger(logDebug, "Selecting elements randomly...")
	std::size_t n = std::distance(begin, end);
	if (k > n) {
		throw std::out_of_range(
				"The number of elements to select is larger than the total number of elements.");
	}
	RandomGen rng(seed);
	for (std::size_t i = 0; i < k; i++) {
		std::size_t ind = rng.getSInt(i, n - 1);
		RandomIterator it = begin + ind;
		std::swap(*(begin + i), *it);
	}
	logger(logDebug, "Done")
}

/**
 * @brief Class for comparing pairs based on second elements only.
 * @tparam FirstT The first pair type.
 * @tparam SecondT The second pair type.
 * @tparam CompareT The comparator type.
 */
template<typename FirstT, typename SecondT, typename CompareT> struct PairCompRight {
	/**
	 * Compare pair based on second elements.
	 * @param l The first pair.
	 * @param r The second pair.
	 * @return The result of the comparison.
	 */
	bool operator()(std::pair<FirstT, SecondT> const & l,
			std::pair<FirstT, SecondT> const & r) const {
		CompareT compare;
		return compare(l.second, r.second);
	}
};

/**
 * Randomly select k largest or smallest elements.
 * @param begin Iterator pointing to the start of elements.
 * @param end Iterator pointing to the end of elements.
 * @param inserter An inserter iterator where the selected elements will be inserted.
 * @param k The number of elements to select. If k is larger than the actual number of elements, std::out_of_range is thrown.
 */
template<typename T, typename InputIterator, typename InserterIt,
		typename Compare> void selectTopK(InputIterator begin,
		InputIterator end, InserterIt inserter, std::size_t k) {
	std::size_t n = std::distance(begin, end);
	if (k > n) {
		throw std::out_of_range(
				"The number of elements to select is larger than the total number of elements.");
	}
	std::priority_queue<T, std::vector<T>, Compare> pq;
	auto it = begin;
	for (std::size_t i = 0; i < k; i++) {
		T val = *it;
		pq.push(val);
		++it;
	}
	Compare compare;
	for (; it != end; ++it) {
		if (!compare(pq.top(), *it)) {
			pq.pop();
			pq.push(*it);
		}
	}
	while (!pq.empty()) {
		inserter = pq.top();
		++inserter;
		pq.pop();
	}
}

/**
 * Select a random element from a range.
 */
template<typename RandomIterator> RandomIterator getRandom(RandomIterator begin,
		RandomIterator end, long int seed) {
	RandomGen rng(seed);
	return begin + rng.getUInt(0, end - begin - 1);
}

/**
 * Filter a range.
 */
template<typename T, typename InputIterator, typename InserterIt> void filter(
		InputIterator begin, InputIterator end, std::set<T> const & excepts,
		InserterIt inserter) {
	for (auto it = begin; it != end; ++it) {
		if (excepts.count(*it) == 0) {
			*inserter = *it;
			++inserter;
		}
	}
}

/**
 * Print data.
 */
template<typename InputIterator> void print(InputIterator begin,
		InputIterator end, std::string const & title, std::ostream & out) {

	out << title << std::endl;
	for (auto it = begin; it != end; ++it) {
		out << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
				<< *it << std::endl;
	}
}

/**
 * Print data.
 */
template<typename InputIterator> void print(InputIterator begin,
		InputIterator end, std::string const & title) {

	print(begin, end, title, std::cout);
}

/**
 * Print edges.
 */
template<typename InputIterator, typename Network> void printEdges(
		InputIterator begin, InputIterator end, Network const & net,
		std::string const & title, std::ostream & out) {

	out << title << std::endl;
	for (auto it = begin; it != end; ++it) {
		out << net.start(*it) << "\t" << net.end(*it) << std::endl;
	}
}

/**
 * Print edges.
 */
template<typename InputIterator, typename Network> void printEdges(
		InputIterator begin, InputIterator end, Network const & net,
		std::string const & title) {

	printEdges(begin, end, net, title, std::cout);
}

/**
 * @return The norm of a range.
 */
template<typename InputIterator> double norm(InputIterator begin,
		InputIterator end) {
	double sum = 0;
	for (auto it = begin; it != end; ++it) {
		sum += (*it) * (*it);
	}
	return std::sqrt(sum);
}

/**
 * Controlled numerical cast.
 */
inline int int_cast(std::size_t source) {
	if (source > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
		throw std::out_of_range(
				"Cannot down cast std::size_t to int, value of the source is out of the range of the target type");
	}
	return static_cast<int>(source);
}

/**
 * Check for NaN and throws exception if it finds it.
 */
// TODO add omp parallel for
template<typename InputIterator> void assertNoNaN(InputIterator begin,
		InputIterator end) {
	for (auto it = begin; it != end; ++it) {
		if (std::isnan(*it)) {
			throw std::runtime_error("Data contains NaN");
		}
	}
}

/**
 * Compute local range in a distributed setting.
 * @param n Total range size.
 * @param nbProcs Number of processors.
 * @param procID The processor ID.
 * @return Pair containing start and end of the range (end is not included in the range).
 */
std::pair<std::size_t, std::size_t> localRange(std::size_t n, int nbProcs,
		int procID);

/**
 * Randomly shuffle a range.
 */
template<typename RandomIterator> void shuffle(RandomIterator begin,
		RandomIterator end, long int seed) {
	std::mt19937_64 g(seed);
	std::shuffle(begin, end, g);
}

#ifdef LINKPRED_WITH_MPI
/**
 * Merge a map-queue over all processes.
 * @param mq Local queue. After call, the queue at process 0 contains the merged queue.
 * @param comm The MPI communicator.
 */
template<typename KComparator, typename PComparator> void merge(
		LMapQueue<unsigned long long, double, KComparator, PComparator> & mq,
		MPI_Comm const & comm) {
	enum Tags {
		MergeTag = 125478
	};

	int nbProcs;
	int procID;
	MPI_Comm_size(comm, &nbProcs);
	if (nbProcs == 1) {
		return; // Nothing to do here
	}
	MPI_Comm_rank(comm, &procID);

// Create an MPI structure for transferring queue content
	int count = 2;
	int blockLengths[2] = {1, 1};
	MPI_Datatype types[2] = {MPI_UNSIGNED_LONG_LONG, MPI_DOUBLE};
	MPI_Datatype mpiPairULLD;
	typedef std::pair<unsigned long long, double> pairULLD;
	MPI_Aint offsets[] =
	{	offsetof(pairULLD, first), offsetof(pairULLD, second)};
	MPI_Type_create_struct(count, blockLengths, offsets, types, &mpiPairULLD);
	MPI_Type_commit(&mpiPairULLD);

	for (int j = 0; j < std::ceil(std::log2(nbProcs)); j++) {
		int d = std::pow(2, j);
		int nd = std::pow(2, j + 1);
		if (procID % nd == 0) { // Receiver
			if (procID + d < nbProcs) {
				std::vector<std::pair<std::size_t, double>> rmq;
				rmq.resize(mq.getL()); // Space for receiving data
				MPI_Status status;
				MPI_Recv(rmq.data(), rmq.size(), mpiPairULLD, procID + d,
						MergeTag, comm, &status);
				int nbRec = 0;
				MPI_Get_count(&status, mpiPairULLD, &nbRec);
				// Merge received data with local queue
				for (int i = nbRec - 1; i >= 0; i--) {
					if (mq.size() == mq.getL()
							&& !mq.compareTop(rmq[i].second)) {
						break;
					}
					mq.push(rmq[i]);
				}
			}
		} else if (procID % d == 0) { // Sender
			// Senders must sort their data
			// Notice that the local queue will be destroyed, but
			// this is OK, senders will not use their local data any more
			std::vector<std::pair<std::size_t, double>> smq;
			if (procID != 0) {
				smq.reserve(mq.size());
				while (!mq.empty()) {
					smq.push_back(mq.top());
					mq.pop();
				}
			}
//			for(auto it = smq.begin(); it != smq.end(); ++it) {
//				std::cout << "Merge: " << it->second << std::endl;
//			}

			MPI_Send(smq.data(), smq.size(), mpiPairULLD, procID - d, MergeTag,
					comm);
		}
		// Wait for all processes to finish
		MPI_Barrier(comm);
	}
}
#endif

/**
 * Read edges from file.
 * @param fileName The file containing the edges.
 * @param ignoreLoops Whether to ignore loops (an edge from a node to itself).
 * @return A vector containing the edge represented as std::pair.
 */
template<typename Label = std::string> std::vector<std::pair<Label, Label>> readEdges(
		std::string fileName, bool ignoreLoops = true) {

	logger(logDebug, "Reading edges from file...")
	std::ifstream in;
	in.open(fileName.c_str(), std::fstream::in);

	if (!in) {
		throw std::runtime_error("Cannot open file: " + fileName);
	}

	std::vector<std::pair<Label, Label>> edges;

	Label start, end;
	std::string line;

	while (std::getline(in, line)) {
		Utils::ltrim(line); // Trim from left.
		if (line.length() == 0 || line[0] == '#') {
			continue;
		}
		std::stringstream sline(line);
		sline >> start;
		sline >> end;
		if ((start != end) || !ignoreLoops) {
			edges.push_back(std::make_pair(start, end));
		}
	}
	in.close();
	logger(logDebug, "Done")
	return edges;
}

/**
 * @brief A structure to store the score of an edge.
 * @tparam Label The label type.
 */
template<typename Label = std::string> struct EdgeScore {
	Label i; /**< The start node. */
	Label j; /**< The end node. */
	double score; /**< The score. */
};

/**
 * Read edges and their scores from file.
 * @param fileName The file containing the edges and their scores.
 * @return A vector containing the edges with their scores.
 */
template<typename Label = std::string> std::vector<EdgeScore<Label>> readEdgeScores(
		std::string fileName) {

	logger(logDebug, "Reading edge scores from file...")
	std::ifstream in;
	in.open(fileName.c_str(), std::fstream::in);

	if (!in) {
		throw std::runtime_error("Cannot open file: " + fileName);
	}

	std::vector<EdgeScore<Label>> edgeScores;

	Label start, end;
	double score;
	std::string line;

	while (std::getline(in, line)) {
		Utils::ltrim(line); // Trim from left.
		if (line.length() == 0 || line[0] == '#') {
			continue;
		}
		std::stringstream sline(line);
		sline >> start;
		sline >> end;
		sline >> score;
		edgeScores.push_back( { start, end, score });
	}
	in.close();
	logger(logDebug, "Done")
	return edgeScores;
}

/**
 * Write edges and their scores to file.
 * @param fileName The file where the edges and their scores will be written.
 * @param esv A vector containing the edges with their scores.
 */
template<typename Label = std::string> void writeEdgeScores(
		std::string fileName, std::vector<EdgeScore<Label>> const & esv) {

	logger(logDebug, "Writing edge scores to file...")
	std::ofstream out;
	out.open(fileName.c_str(), std::fstream::out);

	if (!out) {
		throw std::runtime_error("Cannot open file: " + fileName);
	}

	for (auto it = esv.begin(); it != esv.end(); ++it) {
		out << it->i << "\t" << it->j << "\t" << it->score << std::endl;
	}

	out.close();
	logger(logDebug, "Done")
}

}
/* namepsace Utils */

}
/* namespace LinkPred */

#endif /* INCLUDE_UTILITIES_HPP_ */
