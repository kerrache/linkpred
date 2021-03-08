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

#include <linkpred.hpp>
#include <iostream>
using namespace LinkPred;

int main(int argc, char *argv[]) {

	Vec x = { 1, 0, 1 };
	Vec y = { -1, 1, 2 };

	CosineSim csm;
	std::cout << "CosineSim: " << csm.sim(x, y) << std::endl;

	DotProd dp;
	std::cout << "DotProd: " << dp.sim(x, y) << std::endl;

	L1Sim l1;
	std::cout << "L1Sim: " << l1.sim(x, y) << std::endl;

	L2Sim l2;
	std::cout << "L2Sim: " << l2.sim(x, y) << std::endl;

	LPSim l3(3);
	std::cout << "L3Sim: " << l3.sim(x, y) << std::endl;

	Pearson prs;
	std::cout << "Pearson: " << prs.sim(x, y) << std::endl;

	return 0;
}
