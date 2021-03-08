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

// Original copyright notice
/*
 * fitHRG - fits a hierarchical random graph (hrg) model to data
 * Copyright (C) 2005-2012 Aaron Clauset
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * See http: *www.gnu.org/licenses/gpl.txt for more details.
 *
 * ****************************************************************************************************
 * Author       : Aaron Clauset  ( aaronc@santafe.edu | http: *www.santafe.edu/~aaronc/ )
 * Collaborators: Cristopher Moore and Mark Newman
 * Project      : Hierarchical Random Graphs
 * Location     : University of New Mexico, Dept. of Computer Science AND Santa Fe Institute
 * Created      : 21 June 2006
 * Modified     : 16 June 2007
 *			 : 26 December 2007 (cleaned up for public consumption)
 *			 : 27 May 2012      (workaround for quicksort bug)
 *
 * ****************************************************************************************************
 *
 * This program runs the MCMC with HRG model and the input graph G, samples the likelihoods of edges
 * that are not present in G, and writes out a list of them rank-ordered by their average likelihood
 * under the sampled HRG models. The program includes a convergence criterion that attempts to detect
 * when the MCMC has converged on the equilibrium set.
 *
 * ****************************************************************************************************
 *
 *  See http: *www.santafe.edu/~aaronc/randomgraphs/ for more information, updates to the code, etc.
 */

#include "linkpred/predictors/undirected/uhrgpredictor/hrgds.h"

list::list() {
	x = -1;
	next = nullptr;
}
list::~list() {
}

