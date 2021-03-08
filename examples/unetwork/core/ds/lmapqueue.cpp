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

#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;

int main() {
	LMapQueue<int, double> mq1(5);
	mq1.push(1, 2.1);
	mq1.push(2, 1.5);
	mq1.push(3, 3.2);
	mq1.push(4, 2.7);
	mq1.push(5, 1.9);
	mq1.push(6, 2.0);
	mq1.push(7, 3.4);
	mq1.push(8, 2.1);
	mq1.push(9, 0.8);
	mq1.push(10, 2.9);
	mq1.push(0, 4.2);
//	mq1.printPQ();

	LMapQueue<int, double> mq2(5);
	mq2.push(-1, 2.1);
	mq2.push(-2, 1.5);
	mq2.push(-3, 3.2);
	mq2.push(-4, 2.7);
	mq2.push(-5, 1.9);
	mq2.push(-6, 2.0);
	mq2.push(-7, 3.4);
	mq2.push(-8, 2.1);
	mq2.push(-9, 0.8);
	mq2.push(-10, 2.9);
	mq2.push(0, 4.2);
//	mq2.printPQ();

	std::vector<LMapQueue<int, double>> mqa;
	mqa.push_back(mq1);
	mqa.push_back(mq2);

	auto mq = LMapQueue<int, double>::merge(5, mqa.begin(), mqa.end());
	for (auto it = mq.begin(); it != mq.end(); ++it) {
		std::cout << it->first << "\t" << it->second << std::endl;
	}
	return 0;
}
