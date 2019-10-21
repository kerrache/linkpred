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

#include "linkpred/predictors/uhrgpredictor/dendro_pr.hpp"

interns::interns(const int n, long int seed) :
		seed(seed), rng(seed) {
	q = n;
	count = 0;
	edgelist = new ipair[q];
	indexLUT = new int*[q + 1];
	for (int i = 0; i < (q + 1); i++) {
		indexLUT[i] = new int[2];
		indexLUT[i][0] = indexLUT[i][1] = -1;
	}
}
interns::~interns() {
	delete[] edgelist;
	for (int i = 0; i < (q + 1); i++) {
		delete[] indexLUT[i];
	}
	delete[] indexLUT;
}

// ********************************************************************************************************

// NOTE: Returns an address to another object -- do not deallocate
ipair* interns::getEdge(const int i) {
	return &edgelist[i];
}

// ********************************************************************************************************

// NOTE: Returns an address to another object -- do not deallocate
ipair* interns::getRandomEdge() {
	return &edgelist[(int) (std::floor((double) (q) * rng.getDouble(0, 1)))];
}

// ********************************************************************************************************

bool interns::addEdge(const int new_x, const int new_y,
		const short int new_type) {
	// This function adds a new edge (i,j,t,sp) to the list of internal edges. After checking that the inputs
	// fall in the appropriate range of values, it records the new edgelist index in the indexLUT and then
	// puts the input values into that edgelist location.
	if (count < q && new_x >= 0 && new_x < (q + 1) && new_y >= 0
			&& new_y < (q + 2) && (new_type == LEFT || new_type == RIGHT)) {
		if (new_type == LEFT) {
			indexLUT[new_x][0] = count;
		} else {
			indexLUT[new_x][1] = count;
		}
		edgelist[count].x = new_x;
		edgelist[count].y = new_y;
		edgelist[count].t = new_type;
		count++;
		return true;
	} else {
		return false;
	}
}

// ********************************************************************************************************

void interns::printEdgeList() {
	for (int i = 0; i < q; i++) {
		std::cout << "(" << edgelist[i].x << " " << edgelist[i].y << " ";
		if (edgelist[i].t == LEFT) {
			std::cout << "L) ";
		} else if (edgelist[i].t == RIGHT) {
			std::cout << "R) ";
		} else {
			std::cout << "?) ";
		}
	}
	std::cout << std::endl;
	return;
}

// ********************************************************************************************************

bool interns::swapEdges(const int one_x, const int one_y,
		const short int one_type, const int two_x, const int two_y,
		const short int two_type) {
	// The moves on the dendrogram always swap edges, either of which (or both, or neither) can by internal
	// edges. So, this function mirrors that operation for the internal edgelist && indexLUT.

	int index, jndex, temp;
	bool one_isInternal = false;
	bool two_isInternal = false;

	if (one_x >= 0 && one_x < (q + 1) && two_x >= 0 && two_x < (q + 1)
			&& (two_type == LEFT || two_type == RIGHT) && one_y >= 0
			&& one_y < (q + 2) && two_y >= 0 && two_y < (q + 2)
			&& (one_type == LEFT || one_type == RIGHT)) {

		if (one_type == LEFT) {
			temp = 0;
		} else {
			temp = 1;
		}
		if (indexLUT[one_x][temp] > -1) {
			one_isInternal = true;
		}
		if (two_type == LEFT) {
			temp = 0;
		} else {
			temp = 1;
		}
		if (indexLUT[two_x][temp] > -1) {
			two_isInternal = true;
		}

		if (one_isInternal && two_isInternal) {
			if (one_type == LEFT) {
				index = indexLUT[one_x][0];
			} else {
				index = indexLUT[one_x][1];
			}
			if (two_type == LEFT) {
				jndex = indexLUT[two_x][0];
			} else {
				jndex = indexLUT[two_x][1];
			}
			temp = edgelist[index].y;
			edgelist[index].y = edgelist[jndex].y;
			edgelist[jndex].y = temp;

		} else if (one_isInternal) {
			if (one_type == LEFT) {
				index = indexLUT[one_x][0];
				indexLUT[one_x][0] = -1;
			} else {
				index = indexLUT[one_x][1];
				indexLUT[one_x][1] = -1;
			}
			edgelist[index].x = two_x;
			edgelist[index].t = two_type;
			if (two_type == LEFT) {
				indexLUT[two_x][0] = index;
			} else {
				indexLUT[two_x][1] = index;
			} // add new

		} else if (two_isInternal) {
			if (two_type == LEFT) {
				index = indexLUT[two_x][0];
				indexLUT[two_x][0] = -1;
			} else {
				index = indexLUT[two_x][1];
				indexLUT[two_x][1] = -1;
			}
			edgelist[index].x = one_x;
			edgelist[index].t = one_type;
			if (one_type == LEFT) {
				indexLUT[one_x][0] = index;
			} else {
				indexLUT[one_x][1] = index;
			} // add new
		} else {
		} // else neither is internal

		return true;
	} else {
		return false;
	}
}

elementd::elementd() {
	type = DENDRO;
	logL = p = 0.0;
	e = n = 0;
	label = index = -1;
	M = L = R = NULL;
}
elementd::~elementd() {
}

dendro::~dendro() {
	list *curr, *prev;
	if (g != nullptr) {
		delete g;
		g = nullptr;
	}    // O(m)
	if (internal != nullptr) {
		delete[] internal;
		internal = nullptr;
	}    // O(n)
	if (leaf != nullptr) {
		delete[] leaf;
		leaf = nullptr;
	}    // O(n)
	if (d != nullptr) {
		delete d;
		d = nullptr;
	}    // O(n)
	if (paths != nullptr) {
		for (int i = 0; i < n; i++) {
			curr = paths[i];
			while (curr != nullptr) {
				prev = curr;
				curr = curr->next;
				delete prev;
				prev = nullptr;
			}
			paths[i] = nullptr;
		}
		delete[] paths;
	}
	paths = nullptr;
}

// ********************************************************************************************************

void dendro::binarySearchInsert(elementd* x, elementd* y) {
	if (y->p < x->p) {		// go to left subtree
		if (x->L == nullptr) { // check if left subtree is empty
			x->L = y;		// make x left child
			y->M = x;		// make y parent of child
			return;
		} else {
			binarySearchInsert(x->L, y);
		}
	} else {				// go to right subtree
		if (x->R == nullptr) { // check if right subtree is empty
			x->R = y;		// make x right child
			y->M = x;		// make y parent of child
			return;
		} else {
			binarySearchInsert(x->R, y);
		}
	}
	return;
}

// ********************************************************************************************************

list* dendro::binarySearchFind(const double v) {
	list *head = nullptr, *tail = nullptr, *newlist = nullptr;
	elementd *current = root;
	bool flag_stopSearch = false;

	while (!flag_stopSearch) {				// continue until we're finished
		newlist = new list;				// add this node to the path
		newlist->x = current->label;
		if (current == root) {
			head = newlist;
			tail = head;
		} else {
			tail->next = newlist;
			tail = newlist;
		}
		if (v < current->p) {				// now try left subtree
			if (current->L->type == GRAPH) {
				flag_stopSearch = true;
			} else {
				current = current->L;
			}
		} else {							// else try right subtree
			if (current->R->type == GRAPH) {
				flag_stopSearch = true;
			} else {
				current = current->R;
			}
		}
	}
	return head;
}

// ********************************************************************************************************

void dendro::buildDendrogram() {
	if (g == nullptr) {
		std::cout
				<< "Error: cannot build dendrogram without a graph structure.\n";
		return;
	}

	/* the initialization of the dendrogram structure goes like this:
	 * 1) we allocate space for the n-1 internal nodes of the dendrogram, and then the n leaf nodes
	 * 2) we build a random binary tree structure out of the internal nodes by assigning each
	 *    a uniformly random value over [0,1] and then inserting it into the tree according to the
	 *    binary-search rule.
	 * 3) next, we make a random permutation of the n leaf nodes and add them to the dendrogram D by
	 *    replastd::cing the emptpy spots in-order
	 * 4) then, we compute the path from the root to each leaf and store that in each leaf (this is
	 *    prep work for the next step)
	 * 5) finally, we compute the values for nL, nR, e (and thus p) and the label for each internal
	 *    node by allocating each of the m edges in g to the appropriate internal node
	 */

	// --- Initialization and memory allocation for data structures
	// After allocating the memory for D and G, we need to mark the nodes for G as being
	// non-internal vertices, and then insert them into a random binary tree structure.
	// For simplicity, we make the first internal node in the array the root.
	bool flag_debug = false;
	n = g->getNumNodes();		// size of graph
	leaf = new elementd[n];		// allocate memory for G, O(n)
	internal = new elementd[n - 1];	// allocate memory for D, O(n)
	d = new interns(n - 2, rng.getInt()); // allocate memory for internal edges of D, O(n)
	for (int i = 0; i < n; i++) {		// initialize leaf nodes
		leaf[i].type = GRAPH;
		leaf[i].label = i;
		leaf[i].index = i;
		leaf[i].n = 1;
	}
	if (flag_debug) {
		std::cout << ">> dendro: allocated memory for internal and leaf arrays"
				<< std::endl;
	}
	root = &internal[0];		// initialize internal nodes
	root->label = 0;
	root->index = 0;
	root->p = rng.getDouble(0, 1);
	for (int i = 1; i < (n - 1); i++) {	// insert remaining internal vertices, O(n log n)
		internal[i].label = i;
		internal[i].index = i;
		internal[i].p = rng.getDouble(0, 1);
		binarySearchInsert(root, &internal[i]);
	}
	if (flag_debug) {
		std::cout
				<< ">> dendro: inserted internal vertices into random binary tree"
				<< std::endl;
	}

	// --- Hang leaf nodes off end of dendrogram O(n log n)
	// To impose this random hierarchical relationship on G, we first take a random permutation
	// of the leaf vertices and then replace the nullptrs at the bottom of the tree in-order with
	// the leafs. As a hack to ensure that we can find the leafs later using a binary search,
	// we assign each of them the p value of their parent, perturbed slightly so as to preserve
	// the binary search property.

	block* array;
	array = new block[n];
	for (int i = 0; i < n; i++) {
		array[i].x = rng.getDouble(0, 1);
		array[i].y = i;
	}
	QsortMain(array, 0, n - 1);

	int k = 0;						// replace nullptrs with leaf nodes, and
	for (int i = 0; i < (n - 1); i++) {	//    maintain binary search property, O(n)
		if (internal[i].L == nullptr) {
			internal[i].L = &leaf[array[k].y];
			leaf[array[k].y].M = &internal[i];
			leaf[array[k++].y].p = internal[i].p - 0.0000000000001;
		}
		if (internal[i].R == nullptr) {
			internal[i].R = &leaf[array[k].y];
			leaf[array[k].y].M = &internal[i];
			leaf[array[k++].y].p = internal[i].p + 0.0000000000001;
		}
	}
	delete[] array;
	if (flag_debug) {
		std::cout << ">> dendro: replaced nullptrs in bin-tree with leaf nodes"
				<< std::endl;
	}

	// --- Compute the path from root -> leaf for each leaf O(n log n)
	// Using the binary search property, we can find each leaf node in O(log n) time. The
	// binarySearchFind() function returns the list of internal node indices that the search
	// crossed, in the order of root -> ... -> leaf, for use in the subsequent few operations.

	if (paths != nullptr) {
		list *curr, *prev;
		for (int i = 0; i < n; i++) {
			curr = paths[i];
			while (curr != nullptr) {
				prev = curr;
				curr = curr->next;
				delete prev;
				prev = nullptr;
			}
			paths[i] = nullptr;
		}
		delete[] paths;
	}
	paths = nullptr;
	paths = new list*[n];
	for (int i = 0; i < n; i++) {
		paths[i] = binarySearchFind(leaf[i].p);
	}

	if (flag_debug) {
		std::cout << ">> dendro: computed paths from root to leafs"
				<< std::endl;
	}

	// --- Count e for each internal node O(m)
	// To count the number of edges that span the L and R subtrees for each internal node,
	// we use the path information we just computed. Then, we loop over all edges in G
	// and find the common ancestor in D of the two endpoints and increment that internal
	// node's e count. This process takes O(m) time because in a roughly balanced binary
	// tree (given by our random dendrogram), the vast majority of vertices take basically
	// constant time to find their common ancestor. Note that because our adjacency list
	// is symmetric, we overcount each e by a factor of 2, so we need to correct this after.

	elementd* ancestor;
	edge* curr;
//	int degree;
	for (int i = 0; i < (n - 1); i++) {
		internal[i].e = 0;
		internal[i].label = -1;
	}
	for (int i = 0; i < n; i++) {
		curr = g->getNeighborList(i);
		while (curr != nullptr) {
			ancestor = findCommonAncestor(paths, i, curr->x);
			ancestor->e += 1;
			curr = curr->next;
		}
	}
	for (int i = 0; i < (n - 1); i++) {
		internal[i].e /= 2;
	}
	if (flag_debug) {
		std::cout << ">> dendro: finished common ancestor computation"
				<< std::endl;
	}

	// --- Count n for each internal node O(n log n)
	// To tabulate the number of leafs in each subtree rooted at an internal node,
	// we use the path information computed above.

	for (int i = 0; i < n; i++) {
		ancestor = &leaf[i];
		ancestor = ancestor->M;
		while (ancestor != nullptr) {
			ancestor->n++;
			ancestor = ancestor->M;
		}
	}
	if (flag_debug) {
		std::cout << ">> dendro: computed subtree sizes" << std::endl;
	}

	// --- Label all internal vertices O(n log n)
	// We want to label each internal vertex with the smallest leaf index of its children.
	// This will allow us to collapse many leaf-orderings into a single dendrogram structure
	// that is independent of child-exhanges (since these have no impact on the likelihood
	// of the hierarchical structure). To do this, we loop over the leaf vertices from
	// smallest to largest and walk along that leaf's path from the root. If we find an
	// unlabeled internal node, then we mark it with this leaf's index.

	for (int i = 0; i < n; i++) {
		ancestor = &leaf[i];
		while (ancestor != nullptr) {
			if (ancestor->label == -1 || ancestor->label > leaf[i].label) {
				ancestor->label = leaf[i].label;
			}
			ancestor = ancestor->M;
		}
	}
	if (flag_debug) {
		std::cout << ">> dendro: labeled all internal vertices" << std::endl;
	}

	// --- Exchange children to enforce order-property O(n)
	// We state that the order-property requires that an internal node's label is the
	// smallest index of its left subtree. The dendrogram so far doesn't reflect this, so we
	// need to step through each internal vertex and make that adjustment (swapping nL and nR
	// if we make a change).

//	int temp;
	elementd *tempe;
	for (int i = 0; i < (n - 1); i++) {
		if (internal[i].L->label > internal[i].label) {
			tempe = internal[i].L;
			internal[i].L = internal[i].R;
			internal[i].R = tempe;
		}
	}
	if (flag_debug) {
		std::cout << ">> dendro: enforced order-property" << std::endl;
	}

	// --- Tabulate internal dendrogram edges O(n^2)
	// For the MCMC moves later on, we'll need to be able to choose, uniformly at random, an
	// internal edge of the dendrogram to manipulate. There are always n-2 of them, and we can
	// find them simply by scanning across the internal vertices and observing which have children
	// that are also internal vertices. Note: very important that the order property be enforced
	// before this step is taken; otherwise, the internal edges wont reflect the actual dendrogram
	// structure.

	for (int i = 0; i < (n - 1); i++) {
		if (internal[i].L->type == DENDRO) {
			d->addEdge(i, internal[i].L->index, LEFT);
		}
		if (internal[i].R->type == DENDRO) {
			d->addEdge(i, internal[i].R->index, RIGHT);
		}
	}
	if (flag_debug) {
		std::cout << ">> dendro: tabulated internal dendrogram edges"
				<< std::endl;
	}

	// --- Clear memory for paths O(n log n)
	// Now that we're finished using the paths, we need to deallocate them manually.

	list *current, *previous;
	for (int i = 0; i < n; i++) {
		current = paths[i];
		while (current != nullptr) {
			previous = current;
			current = current->next;
			delete previous;
			previous = nullptr;
		}
		paths[i] = nullptr;
	}
	delete[] paths;
	paths = nullptr;
	if (flag_debug) {
		std::cout << ">> dendro: cleared memory for paths" << std::endl;
	}

	// --- Compute p_i for each internal node O(n)
	// Each internal node's p_i = e_i / (nL_i*nR_i), and now that we have each of those
	// pieces, we may calculate this value for each internal node. Given these, we can then
	// calculate the log-likelihood of the entire dendrogram structure
	// \log(L) = \sum_{i=1}^{n} ( ( e_i \log[p_i] ) + ( (nL_i*nR_i - e_i) \log[1-p_i] ) )

	L = 0.0;
	double dL;
	int nL_nR, ei;
	for (int i = 0; i < (n - 1); i++) {
		nL_nR = internal[i].L->n * internal[i].R->n;
		ei = internal[i].e;
		internal[i].p = (double) (ei) / (double) (nL_nR);
		if (ei == 0 || ei == nL_nR) {
			dL = 0.0;
		} else {
			dL = ei * log(internal[i].p)
					+ (nL_nR - ei) * log(1.0 - internal[i].p);
		}
		internal[i].logL = dL;
		L += dL;
	}

	char pauseme;
	for (int i = 0; i < (n - 1); i++) {
		if (internal[i].label > internal[i].L->label) {
			tempe = internal[i].L;
			internal[i].L = internal[i].R;
			internal[i].R = tempe;
			std::cout << "#### WARNING - order property violated by internal["
					<< i << "] (fixed)" << std::endl;
			std::cin >> pauseme;
		}
	}

	// --- Dendrogram is now built
	if (flag_debug) {
		std::cout << ">> dendro: build dendrogram complete: logL = " << L
				<< std::endl;
	}

	return;
}

// ********************************************************************************************************

void dendro::clearDendrograph() {
	// Clear out the memory and references used by the dendrograph structure - this is
	// intended to be called just before an importDendrogramStructure call so as to avoid
	// memory leaks and overwriting the references therein.
	if (g != nullptr) {
		delete g;
		g = nullptr;
	}    // O(m)
	if (leaf != nullptr) {
		delete[] leaf;
		leaf = nullptr;
	}    // O(n)
	if (internal != nullptr) {
		delete[] internal;
		internal = nullptr;
	}    // O(n)
	if (d != nullptr) {
		delete d;
		d = nullptr;
	}    // O(n)
	root = nullptr;

	return;
}

// ********************************************************************************************************

int dendro::computeEdgeCount(const int a, const short int atype, const int b,
		const short int btype) {
	// This function computes the number of edges that cross between the subtree internal[a]
	// and the subtree internal[b].
	// To do this, we use an array A[1..n] integers which take values -1 if A[i] is in the
	// subtree defined by internal[a], +1 if A[i] is in the subtree internal[b], and 0
	// otherwise. Taking the smaller of the two sets, we then scan over the edges attached
	// to that set of vertices and count the number of endpoints we see in the other set.

//	bool flag_debug = true;
	bool flag_go = true;
	int nA, nB;	// q;
	int count = 0;
	const short int k = 1 + DENDRO + GRAPH;

	elementd* curr;

	// --- First, we push the leaf nodes in the L and R subtrees into balanced binary tree
	//     structures so that we can search them quickly later on.
	if (atype == GRAPH) {					// default case, subtree A is size 1
		subtreeL.insertItem(a, -1);	// insert single node as member of left subtree
		nA = 1;					//
	} else {
		curr = &internal[a];			// explore subtree A, O(|A|)
		curr->type = k + 1;					//
		nA = 0;					//
		while (flag_go) {

			if (curr->index == internal[a].M->index) {
				internal[a].type = DENDRO;
				flag_go = false;
			} else {
				if (curr->type == k + 1 && curr->L->type == GRAPH) {// - is it time, and is left child a graph node?
					subtreeL.insertItem(curr->L->index, -1);
					curr->type = k + 2;  //
					nA++;					//
				}
				if (curr->type == k + 2 && curr->R->type == GRAPH) {// - is it time, and is right child a graph node?
					subtreeL.insertItem(curr->R->index, -1);
					curr->type = k + 3;  //
					nA++;					//
				}
				if (curr->type == k + 1) {			// - go left
					curr->type = k + 2;			//
					curr = curr->L;		//
					curr->type = k + 1;
				} else if (curr->type == k + 2) {		// - else go right
					curr->type = k + 3;			//
					curr = curr->R;		//
					curr->type = k + 1;
				} else {						// - else go up a level
					curr->type = DENDRO;		//
					curr = curr->M;		//
					if (curr == nullptr) {
						flag_go = false;
					}
				}
			}
			if (nA > n) {
				std::cout
						<< "error! nA > n\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n"
						<< std::endl;
				break;
			}
		}
	}

	if (btype == GRAPH) {					// default case, subtree A is size 1
		subtreeR.insertItem(b, 1);// insert node as single member of right subtree
		nB = 1;					//
	} else {
		flag_go = true;
		curr = &internal[b];			// explore subtree B, O(|B|)
		curr->type = k + 1;					//
		nB = 0;					//
		while (flag_go) {

			if (curr->index == internal[b].M->index) {
				internal[b].type = DENDRO;
				flag_go = false;
			} else {
				if (curr->type == k + 1 && curr->L->type == GRAPH) {// - is it time, and is left child a graph node?
					subtreeR.insertItem(curr->L->index, 1);
					curr->type = k + 2;  //
					nB++;					//
				}
				if (curr->type == k + 2 && curr->R->type == GRAPH) {// - is it time, and is right child a graph node?
					subtreeR.insertItem(curr->R->index, 1);
					curr->type = k + 3;  //
					nB++;					//
				}
				if (curr->type == k + 1) {			// - look left
					curr->type = k + 2;			//
					curr = curr->L;		//
					curr->type = k + 1;
				} else if (curr->type == k + 2) {		// - look right
					curr->type = k + 3;			//
					curr = curr->R;		//
					curr->type = k + 1;
				} else {						// - else go up a level
					curr->type = DENDRO;		//
					curr = curr->M;		//
					if (curr == nullptr) {
						flag_go = false;
					}
				}
			}
			if (nB > n) {
				std::cout
						<< "error! nB > n \n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n#\n"
						<< std::endl;
				break;
			}
		}
	}

	// --- Now, we take the smaller subtree and ask how many of its emerging edges have their
	//     partner in the other subtree. O(|A| log |A|) time
	edge* current;
	int* treeList;
	if (nA < nB) {
		treeList = subtreeL.returnArrayOfKeys();	// subtreeL is smaller
		for (int i = 0; i < nA; i++) {
			current = g->getNeighborList(treeList[i]);
			while (current != nullptr) {// loop over each of its neighbors v_j
				if (subtreeR.findItem(current->x) != nullptr) {
					count++;
				}
				current = current->next;			// to see if v_j is in A
			}								//
			subtreeL.deleteItem(treeList[i]);
		}
		delete[] treeList;
		treeList = subtreeR.returnArrayOfKeys();
		for (int i = 0; i < nB; i++) {
			subtreeR.deleteItem(treeList[i]);
		}
		delete[] treeList;
	} else {
		treeList = subtreeR.returnArrayOfKeys();	// subtreeR is smaller
		for (int i = 0; i < nB; i++) {
			current = g->getNeighborList(treeList[i]);
			while (current != nullptr) {// loop over each of its neighbors v_j
				if (subtreeL.findItem(current->x) != nullptr) {
					count++;
				}
				current = current->next;			// to see if v_j is in B
			}								//
			subtreeR.deleteItem(treeList[i]);
		}
		delete[] treeList;
		treeList = subtreeL.returnArrayOfKeys();
		for (int i = 0; i < nA; i++) {
			subtreeL.deleteItem(treeList[i]);
		}
		delete[] treeList;
	}

	return count;
}

// ********************************************************************************************************

elementd* dendro::findCommonAncestor(list** paths, const int i, const int j) {
	list* headOne = paths[i];
	list* headTwo = paths[j];
	elementd* lastStep = nullptr;
	while (headOne->x == headTwo->x) {
		lastStep = &internal[headOne->x];
		headOne = headOne->next;
		headTwo = headTwo->next;
		if (headOne == nullptr || headTwo == nullptr) {
			break;
		}
	}
	return lastStep;// Returns address of an internal node; do not deallocate
}

// ********************************************************************************************************

double dendro::getLikelihood() {
	if (L > 0.0) {
		return 0.0;
	} else {
		return L;
	}
}

// ********************************************************************************************************

bool dendro::importDendrogramStructure(const std::string in_file) {
	std::string bracketL, bracketR, sL, sR, sLtype, sRtype, sp, se, sn;
	int sindex, sLindex, sRindex, snume, snumn;
	double sprob;
	bool safeExit = true;
	bool flag_debug = true;
	n = 1;

	std::cout << ">> dendro: importing hrg file: " << in_file << std::endl;
	std::ifstream fscan(in_file.c_str(), std::ios::in);
	while (fscan >> sn) {
		if (sn == "[") {
			n++;
		}
	}
	fscan.close();

	leaf = new elementd[n];		// allocate memory for G, O(n)
	internal = new elementd[n - 1];	// allocate memory for D, O(n)
	d = new interns(n - 2, rng.getInt());// allocate memory for internal edges of D, O(n)
	for (int i = 0; i < n; i++) {		// initialize leaf nodes
		leaf[i].type = GRAPH;
		leaf[i].label = i;
		leaf[i].index = i;
		leaf[i].n = 1;
	}
	root = &internal[0];		// initialize internal nodes
	root->label = 0;
	for (int i = 1; i < (n - 1); i++) {
		internal[i].index = i;
		internal[i].label = -1;
	}
	if (flag_debug) {
		std::cout << ">> dendro: allocated memory for internal and leaf arrays"
				<< std::endl;
	}

	// --- Import basic structure froms file O(n)
	std::ifstream fin(in_file.c_str(), std::ios::in);
	while (fin >> bracketL >> sindex >> bracketR >> sL >> sLindex >> sLtype
			>> sR >> sRindex >> sRtype >> sp >> sprob >> se >> snume >> sn
			>> snumn) {
		std::cout << bracketL << " " << sindex << " " << bracketR << " " << sL
				<< " " << sLindex << " " << sLtype << " " << sR << " "
				<< sRindex << " " << sRtype << " " << sp << " " << sprob << " "
				<< se << " " << snume << " " << sn << " " << snumn << std::endl;
		if (sLtype == "(D)") {
			internal[sindex].L = &internal[sLindex];
			internal[sLindex].M = &internal[sindex];
		} else if (sLtype == "(G)") {
			internal[sindex].L = &leaf[sLindex];
			leaf[sLindex].M = &internal[sindex];
		} else {
			std::cout << "Error: " << bracketL << sindex << bracketR << sL
					<< sLindex << sLtype << sR << sRindex << sRtype << sp
					<< sprob << se << snume << sn << snumn << std::endl;
			safeExit = false;
			break;
		}
		if (sRtype == "(D)") {
			internal[sindex].R = &internal[sRindex];
			internal[sRindex].M = &internal[sindex];
		} else if (sRtype == "(G)") {
			internal[sindex].R = &leaf[sRindex];
			leaf[sRindex].M = &internal[sindex];
		} else {
			std::cout << "Error: " << bracketL << sindex << bracketR << sL
					<< sLindex << sLtype << sR << sRindex << sRtype << sp
					<< sprob << se << snume << sn << snumn << std::endl;
			safeExit = false;
			break;
		}
		internal[sindex].p = sprob;
		if (sprob < 0.0 || sprob > 1.0) {
			std::cout << "Error: " << bracketL << sindex << bracketR << sL
					<< sLindex << sLtype << sR << sRindex << sRtype << sp
					<< sprob << se << snume << sn << snumn << std::endl;
			safeExit = false;
			break;
		}
		internal[sindex].e = snume;
		internal[sindex].n = snumn;
		internal[sindex].index = sindex;
	}
	fin.close();
	if (!safeExit) {
		return false;
	}
	if (flag_debug) {
		std::cout << ">> dendro: imported basic structure" << std::endl;
	}

	// --- Label all internal vertices O(n log n)
	elementd* curr;
	for (int i = 0; i < n; i++) {
		curr = &leaf[i];
		while (curr != nullptr) {
			if (curr->label == -1 || curr->label > leaf[i].label) {
				curr->label = leaf[i].label;
			}
			curr = curr->M;
		}
	}
	if (flag_debug) {
		std::cout << ">> dendro: labeled all internal vertices" << std::endl;
	}

	// --- Exchange children to enforce order-property O(n)
//	int temp;
	elementd *tempe;
	for (int i = 0; i < (n - 1); i++) {
		if (internal[i].L->label > internal[i].label) {
			tempe = internal[i].L;
			internal[i].L = internal[i].R;
			internal[i].R = tempe;
		}
	}
	if (flag_debug) {
		std::cout << ">> dendro: enforced order-property" << std::endl;
	}

	// --- Tabulate internal dendrogram edges O(n)
//	int k = 0;
	for (int i = 0; i < (n - 1); i++) {
		if (internal[i].L->type == DENDRO) {
			d->addEdge(i, internal[i].L->index, LEFT);
		}
		if (internal[i].R->type == DENDRO) {
			d->addEdge(i, internal[i].R->index, RIGHT);
		}
	}
	if (flag_debug) {
		std::cout << ">> dendro: tabulated internal dendrogram edges"
				<< std::endl;
	}

	// --- Compute p_i for each internal node O(n)
	// Each internal node's p_i = e_i / (nL_i*nR_i), and now that we have each of those
	// pieces, we may calculate this value for each internal node. Given these, we can then
	// calculate the log-likelihood of the entire dendrogram structure
	// \log(L) = \sum_{i=1}^{n} ( ( e_i \log[p_i] ) + ( (nL_i*nR_i - e_i) \log[1-p_i] ) )

	L = 0.0;
	double dL;
	int nL_nR, ei;
	for (int i = 0; i < (n - 1); i++) {
		nL_nR = internal[i].L->n * internal[i].R->n;
		ei = internal[i].e;
		if (ei == 0 || ei == nL_nR) {
			dL = 0.0;
		} else {
			dL = (double) (ei) * log(internal[i].p)
					+ (double) (nL_nR - ei) * log(1.0 - internal[i].p);
		}
		internal[i].logL = dL;
		L += dL;
	}
	if (flag_debug) {
		std::cout << ">> dendro: computed log-likelihood" << std::endl;
		std::cout << "Log-Likelihood = " << L << std::endl;
	}

	// --- Dendrogram is now built
	if (flag_debug) {
		std::cout << ">> dendro: build dendrogram complete" << std::endl;
	}

	return true;
}

// ********************************************************************************************************

bool dendro::monteCarloMove(double& delta, bool& ftaken) {
	// A single MC move begins with the selection of a random internal edge (a,b) of the
	// dendrogram. This also determines the three subtrees i, j, k that we will rearrange,
	// and we choose uniformly from among the options.
	//
	// If (a,b) is a left-edge, then we have ((i,j),k), and moves
	// ((i,j),k) -> ((i,k),j)								(alpha move)
	//           -> (i,(j,k)) + enforce order-property for (j,k)	(beta move)
	//
	// If (a,b) is a right-edge, then we have (i,(j,k)), and moves
	// (i,(j,k)) -> ((i,k),j)								(alpha move)
	//           -> ((i,j),k)								(beta move)
	//
	// For each of these moves, we need to know what the change in likelihood will be, so
	// that we can determine with what probability we execute the move.

	bool flag_debug = true;
	elementd *temp;	// *tempe;
	ipair *tempPair;
	int x, y, e_x, e_y, n_i, n_j, n_k, n_x, n_y;	// tempi;
	short int t;
	double p_x, p_y, L_x, L_y, dLogL;
//	char pauseme;

	// The remainder of the code executes a single MCMC move, where we sample the dendrograms
	// proportionally to their likelihoods (i.e., temperature=1, if you're comparing it to the
	// usual MCMC framework).
	delta = 0.0;
	ftaken = false;
	tempPair = d->getRandomEdge();	// returns address; no need to deallocate
	x = tempPair->x;			// copy contents of referenced random edge
	y = tempPair->y;			//    into local variables
	t = tempPair->t;

	if (flag_debug) {
		if (t == LEFT) {
		} else if (t == RIGHT) {
		} else {
			std::cout << " bad edge (i,j,k)" << std::endl;
		}
	}

	if (t == LEFT) {								//
		if (rng.getDouble(0, 1) < 0.5) {// ## LEFT ALPHA move: ((i,j),k) -> ((i,k),j)
			// We need to calculate the change in the likelihood that would result from
			// this move. Most of the information needed to do this is already available,
			// the exception being e_ik, the number of edges that span the i and k subtrees.
			// We use a slow algorithm O(n) to do this, since I don't know of a better way at
			// this point.

			n_i = internal[y].L->n;
			n_j = internal[y].R->n;
			n_k = internal[x].R->n;

			n_y = n_i * n_k;
			e_y = computeEdgeCount(internal[y].L->index, internal[y].L->type,
					internal[x].R->index, internal[x].R->type);   // e_ik
			p_y = (double) (e_y) / (double) (n_y);
			if (e_y == 0 || e_y == n_y) {
				L_y = 0.0;
			} else {
				L_y = (double) (e_y) * log(p_y)
						+ (double) (n_y - e_y) * log(1.0 - p_y);
			}

			n_x = (n_i + n_k) * n_j;
			e_x = internal[x].e + internal[y].e - e_y;					// e_yj
			p_x = (double) (e_x) / (double) (n_x);
			if (e_x == 0 || e_x == n_x) {
				L_x = 0.0;
			} else {
				L_x = (double) (e_x) * log(p_x)
						+ (double) (n_x - e_x) * log(1.0 - p_x);
			}

			dLogL = (L_x - internal[x].logL) + (L_y - internal[y].logL);
			if ((dLogL > 0.0) || (rng.getDouble(0, 1) < exp(dLogL))) { // make LEFT ALPHA move
				ftaken = true;
				d->swapEdges(x, internal[x].R->index, RIGHT, y,
						internal[y].R->index, RIGHT);
				temp = internal[x].R;			// - swap j and k
				internal[x].R = internal[y].R;			//
				internal[y].R = temp;					//
				internal[x].R->M = &internal[x];	// - adjust parent pointers
				internal[y].R->M = &internal[y];			//
				internal[y].n = n_i + n_k;				// - update n for [y]
				internal[x].e = e_x;			// - update e_i for [x] and [y]
				internal[y].e = e_y;					//
				internal[x].p = p_x;			// - update p_i for [x] and [y]
				internal[y].p = p_y;					//
				internal[x].logL = L_x;			// - update L_i for [x] and [y]
				internal[y].logL = L_y;					//
				// - order-property maintained
				L += dLogL;				// - update LogL
				delta = dLogL;				//

				// TRAP: Catches violations of the ordering property
				if (internal[x].label > internal[x].L->label
						|| internal[y].label > internal[y].L->label) {
					printDendrogram();
					if (internal[x].label > internal[x].L->label) {
						std::cout
								<< "**** WARNING - order property violated by internal["
								<< x << "]" << std::endl;
						std::cout << "x    (p = " << internal[x].p << "\te = "
								<< internal[x].e << "\tnL = ";
						std::cout << internal[x].L->n << "\tnR = "
								<< internal[x].R->n << "\tlabel = "
								<< internal[x].label << ")\tinternal[" << x
								<< "]\t(D)" << std::endl;
						if (internal[x].L->type == GRAPH) {
							std::cout << "x->L [" << internal[x].L->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->L (p = " << internal[x].L->p
									<< "\te = " << internal[x].L->e
									<< "\tnL = ";
							std::cout << internal[x].L->L->n << "\tnR = "
									<< internal[x].L->R->n << "\tlabel = "
									<< internal[x].L->label << ")\tinternal["
									<< internal[x].L->index << "]\t(D)"
									<< std::endl;
						}
						if (internal[x].R->type == GRAPH) {
							std::cout << "x->R [" << internal[x].R->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->R (p = " << internal[x].R->p
									<< "\te = " << internal[x].R->e
									<< "\tnL = ";
							std::cout << internal[x].R->L->n << "\tnR = "
									<< internal[x].R->R->n << "\tlabel = "
									<< internal[x].R->label << ")\tinternal["
									<< internal[x].R->index << "]\t(D)"
									<< std::endl;
						}
					}
					if (internal[y].label > internal[y].L->label) {
						std::cout
								<< "**** WARNING - order property violated by internal["
								<< y << "]" << std::endl;
						std::cout << "y    (p = " << internal[y].p << "\te = "
								<< internal[y].e << "\tnL = ";
						std::cout << internal[y].L->n << "\tnR = "
								<< internal[y].R->n << "\tlabel = "
								<< internal[y].label << ")\tinternal[" << y
								<< "]\t(D)" << std::endl;
						if (internal[y].L->type == GRAPH) {
							std::cout << "y->L [" << internal[y].L->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->L (p = " << internal[y].L->p
									<< "\te = " << internal[y].L->e
									<< "\tnL = ";
							std::cout << internal[y].L->L->n << "\tnR = "
									<< internal[y].L->R->n << "\tlabel = "
									<< internal[y].L->label << ")\tinternal["
									<< internal[y].L->index << "]\t(D)"
									<< std::endl;
						}
						if (internal[y].R->type == GRAPH) {
							std::cout << "y->R [" << internal[y].R->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->R (p = " << internal[y].R->p
									<< "\te = " << internal[y].R->e
									<< "\tnL = ";
							std::cout << internal[y].R->L->n << "\tnR = "
									<< internal[y].R->R->n << "\tlabel = "
									<< internal[y].R->label << ")\tinternal["
									<< internal[y].R->index << "]\t(D)"
									<< std::endl;
						}
					}
					return false;
				}
			}
		} else {				// ## LEFT BETA move:  ((i,j),k) -> (i,(j,k))
			n_i = internal[y].L->n;
			n_j = internal[y].R->n;
			n_k = internal[x].R->n;

			n_y = n_j * n_k;
			e_y = computeEdgeCount(internal[y].R->index, internal[y].R->type,
					internal[x].R->index, internal[x].R->type);   // e_jk
			p_y = (double) (e_y) / (double) (n_y);
			if (e_y == 0 || e_y == n_y) {
				L_y = 0.0;
			} else {
				L_y = (double) (e_y) * log(p_y)
						+ (double) (n_y - e_y) * log(1.0 - p_y);
			}

			n_x = (n_j + n_k) * n_i;
			e_x = internal[x].e + internal[y].e - e_y;					// e_yj
			p_x = (double) (e_x) / (double) (n_x);
			if (e_x == 0 || e_x == n_x) {
				L_x = 0.0;
			} else {
				L_x = (double) (e_x) * log(p_x)
						+ (double) (n_x - e_x) * log(1.0 - p_x);
			}

			dLogL = (L_x - internal[x].logL) + (L_y - internal[y].logL);
			if ((dLogL > 0.0) || (rng.getDouble(0, 1) < exp(dLogL))) { // make LEFT BETA move
				ftaken = true;
				d->swapEdges(y, internal[y].L->index, LEFT, y,
						internal[y].R->index, RIGHT);
				temp = internal[y].L;			// - swap L and R of [y]
				internal[y].L = internal[y].R;			//
				internal[y].R = temp;					//
				d->swapEdges(x, internal[x].R->index, RIGHT, y,
						internal[y].R->index, RIGHT);
				temp = internal[x].R;			// - swap i and k
				internal[x].R = internal[y].R;			//
				internal[y].R = temp;					//
				internal[x].R->M = &internal[x];	// - adjust parent pointers
				internal[y].R->M = &internal[y];			//
				d->swapEdges(x, internal[x].L->index, LEFT, x,
						internal[x].R->index, RIGHT);
				temp = internal[x].L;			// - swap L and R of [x]
				internal[x].L = internal[x].R;			//
				internal[x].R = temp;					//
				internal[y].n = n_j + n_k;				// - update n
				internal[x].e = e_x;					// - update e_i
				internal[y].e = e_y;					//
				internal[x].p = p_x;					// - update p_i
				internal[y].p = p_y;					//
				internal[x].logL = L_x;					// - update logL_i
				internal[y].logL = L_y;					//
				if (internal[y].R->label < internal[y].L->label) {
					d->swapEdges(y, internal[y].L->index, LEFT, y,
							internal[y].R->index, RIGHT);
					temp = internal[y].L; // - enforce order-property if necessary
					internal[y].L = internal[y].R;    //
					internal[y].R = temp;			//
				}									//
				internal[y].label = internal[y].L->label;    //
				L += dLogL;				// - update LogL
				delta = dLogL;				//

				// TRAP: Catches violations of the ordering property
				if (internal[x].label > internal[x].L->label
						|| internal[y].label > internal[y].L->label) {
					printDendrogram();
					if (internal[x].label > internal[x].L->label) {
						std::cout
								<< "**** WARNING - order property violated by internal["
								<< x << "]" << std::endl;
						std::cout << "x    (p = " << internal[x].p << "\te = "
								<< internal[x].e << "\tnL = ";
						std::cout << internal[x].L->n << "\tnR = "
								<< internal[x].R->n << "\tlabel = "
								<< internal[x].label << ")\tinternal[" << x
								<< "]\t(D)" << std::endl;
						if (internal[x].L->type == GRAPH) {
							std::cout << "x->L [" << internal[x].L->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->L (p = " << internal[x].L->p
									<< "\te = " << internal[x].L->e
									<< "\tnL = ";
							std::cout << internal[x].L->L->n << "\tnR = "
									<< internal[x].L->R->n << "\tlabel = "
									<< internal[x].L->label << ")\tinternal["
									<< internal[x].L->index << "]\t(D)"
									<< std::endl;
						}
						if (internal[x].R->type == GRAPH) {
							std::cout << "x->R [" << internal[x].R->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->R (p = " << internal[x].R->p
									<< "\te = " << internal[x].R->e
									<< "\tnL = ";
							std::cout << internal[x].R->L->n << "\tnR = "
									<< internal[x].R->R->n << "\tlabel = "
									<< internal[x].R->label << ")\tinternal["
									<< internal[x].R->index << "]\t(D)"
									<< std::endl;
						}
					}
					if (internal[y].label > internal[y].L->label) {
						std::cout
								<< "**** WARNING - order property violated by internal["
								<< y << "]" << std::endl;
						std::cout << "y    (p = " << internal[y].p << "\te = "
								<< internal[y].e << "\tnL = ";
						std::cout << internal[y].L->n << "\tnR = "
								<< internal[y].R->n << "\tlabel = "
								<< internal[y].label << ")\tinternal[" << y
								<< "]\t(D)" << std::endl;
						if (internal[y].L->type == GRAPH) {
							std::cout << "y->L [" << internal[y].L->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->L (p = " << internal[y].L->p
									<< "\te = " << internal[y].L->e
									<< "\tnL = ";
							std::cout << internal[y].L->L->n << "\tnR = "
									<< internal[y].L->R->n << "\tlabel = "
									<< internal[y].L->label << ")\tinternal["
									<< internal[y].L->index << "]\t(D)"
									<< std::endl;
						}
						if (internal[y].R->type == GRAPH) {
							std::cout << "y->R [" << internal[y].R->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->R (p = " << internal[y].R->p
									<< "\te = " << internal[y].R->e
									<< "\tnL = ";
							std::cout << internal[y].R->L->n << "\tnR = "
									<< internal[y].R->R->n << "\tlabel = "
									<< internal[y].R->label << ")\tinternal["
									<< internal[y].R->index << "]\t(D)"
									<< std::endl;
						}
					}
					return false;
				}
			}
		}
	} else {										// right-edge: t == RIGHT
		if (rng.getDouble(0, 1) < 0.5) {// alpha move: (i,(j,k)) -> ((i,k),j)
			n_i = internal[x].L->n;
			n_j = internal[y].L->n;
			n_k = internal[y].R->n;

			n_y = n_i * n_k;
			e_y = computeEdgeCount(internal[x].L->index, internal[x].L->type,
					internal[y].R->index, internal[y].R->type);   // e_ik
			p_y = (double) (e_y) / (double) (n_y);
			if (e_y == 0 || e_y == n_y) {
				L_y = 0.0;
			} else {
				L_y = (double) (e_y) * log(p_y)
						+ (double) (n_y - e_y) * log(1.0 - p_y);
			}

			n_x = (n_i + n_k) * n_j;
			e_x = internal[x].e + internal[y].e - e_y;					// e_yj
			p_x = (double) (e_x) / (double) (n_x);
			if (e_x == 0 || e_x == n_x) {
				L_x = 0.0;
			} else {
				L_x = (double) (e_x) * log(p_x)
						+ (double) (n_x - e_x) * log(1.0 - p_x);
			}

			dLogL = (L_x - internal[x].logL) + (L_y - internal[y].logL);
			if ((dLogL > 0.0) || (rng.getDouble(0, 1) < exp(dLogL))) { // make RIGHT ALPHA move
				ftaken = true;
				d->swapEdges(x, internal[x].L->index, LEFT, x,
						internal[x].R->index, RIGHT);
				temp = internal[x].L;			// - swap L and R of [x]
				internal[x].L = internal[x].R;			//
				internal[x].R = temp;				//
				d->swapEdges(y, internal[y].L->index, LEFT, x,
						internal[x].R->index, RIGHT);
				temp = internal[y].L;			// - swap i and j
				internal[y].L = internal[x].R;			//
				internal[x].R = temp;				//
				internal[x].R->M = &internal[x];	// - adjust parent pointers
				internal[y].L->M = &internal[y];			//
				internal[y].n = n_i + n_k;			// - update n
				internal[x].e = e_x;					// - update e_i
				internal[y].e = e_y;					//
				internal[x].p = p_x;					// - update p_i
				internal[y].p = p_y;					//
				internal[x].logL = L_x;					// - update logL_i
				internal[y].logL = L_y;					//
				internal[y].label = internal[x].label;// - update order property
				L += dLogL;				// - update LogL
				delta = dLogL;				//

				// TRAP: Catches violations of the ordering property
				if (internal[x].label > internal[x].L->label
						|| internal[y].label > internal[y].L->label) {
					printDendrogram();
					if (internal[x].label > internal[x].L->label) {
						std::cout
								<< "**** WARNING - order property violated by internal["
								<< x << "]" << std::endl;
						std::cout << "x    (p = " << internal[x].p << "\te = "
								<< internal[x].e << "\tnL = ";
						std::cout << internal[x].L->n << "\tnR = "
								<< internal[x].R->n << "\tlabel = "
								<< internal[x].label << ")\tinternal[" << x
								<< "]\t(D)" << std::endl;
						if (internal[x].L->type == GRAPH) {
							std::cout << "x->L [" << internal[x].L->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->L (p = " << internal[x].L->p
									<< "\te = " << internal[x].L->e
									<< "\tnL = ";
							std::cout << internal[x].L->L->n << "\tnR = "
									<< internal[x].L->R->n << "\tlabel = "
									<< internal[x].L->label << ")\tinternal["
									<< internal[x].L->index << "]\t(D)"
									<< std::endl;
						}
						if (internal[x].R->type == GRAPH) {
							std::cout << "x->R [" << internal[x].R->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->R (p = " << internal[x].R->p
									<< "\te = " << internal[x].R->e
									<< "\tnL = ";
							std::cout << internal[x].R->L->n << "\tnR = "
									<< internal[x].R->R->n << "\tlabel = "
									<< internal[x].R->label << ")\tinternal["
									<< internal[x].R->index << "]\t(D)"
									<< std::endl;
						}
					}
					if (internal[y].label > internal[y].L->label) {
						std::cout
								<< "**** WARNING - order property violated by internal["
								<< y << "]" << std::endl;
						std::cout << "y    (p = " << internal[y].p << "\te = "
								<< internal[y].e << "\tnL = "
								<< internal[y].L->n << "\tnR = "
								<< internal[y].R->n << "\tlabel = "
								<< internal[y].label << ")\tinternal[" << y
								<< "]\t(D)" << std::endl;
						if (internal[y].L->type == GRAPH) {
							std::cout << "y->L [" << internal[y].L->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->L (p = " << internal[y].L->p
									<< "\te = " << internal[y].L->e
									<< "\tnL = ";
							std::cout << internal[y].L->L->n << "\tnR = "
									<< internal[y].L->R->n << "\tlabel = "
									<< internal[y].L->label << ")\tinternal["
									<< internal[y].L->index << "]\t(D)"
									<< std::endl;
						}
						if (internal[y].R->type == GRAPH) {
							std::cout << "y->R [" << internal[y].R->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->R (p = " << internal[y].R->p
									<< "\te = " << internal[y].R->e
									<< "\tnL = ";
							std::cout << internal[y].R->L->n << "\tnR = "
									<< internal[y].R->R->n << "\tlabel = "
									<< internal[y].R->label << ")\tinternal["
									<< internal[y].R->index << "]\t(D)"
									<< std::endl;
						}
					}
					return false;
				}
			}
		} else {						// beta move:  (i,(j,k)) -> ((i,j),k)
			n_i = internal[x].L->n;
			n_j = internal[y].L->n;
			n_k = internal[y].R->n;

			n_y = n_i * n_j;
			e_y = computeEdgeCount(internal[x].L->index, internal[x].L->type,
					internal[y].L->index, internal[y].L->type);   // e_ij
			p_y = (double) (e_y) / (double) (n_y);
			if (e_y == 0 || e_y == n_y) {
				L_y = 0.0;
			} else {
				L_y = (double) (e_y) * log(p_y)
						+ (double) (n_y - e_y) * log(1.0 - p_y);
			}

			n_x = (n_i + n_j) * n_k;
			e_x = internal[x].e + internal[y].e - e_y;					// e_yk
			p_x = (double) (e_x) / (double) (n_x);
			if (e_x == 0 || e_x == n_x) {
				L_x = 0.0;
			} else {
				L_x = (double) (e_x) * log(p_x)
						+ (double) (n_x - e_x) * log(1.0 - p_x);
			}

			dLogL = (L_x - internal[x].logL) + (L_y - internal[y].logL);
			if ((dLogL > 0.0) || (rng.getDouble(0, 1) < exp(dLogL))) { // make RIGHT BETA move
				ftaken = true;
				d->swapEdges(x, internal[x].L->index, LEFT, x,
						internal[x].R->index, RIGHT);
				temp = internal[x].L;			// - swap L and R of [x]
				internal[x].L = internal[x].R;			//
				internal[x].R = temp;				//
				d->swapEdges(x, internal[x].R->index, RIGHT, y,
						internal[y].R->index, RIGHT);
				temp = internal[x].R;			// - swap i and k
				internal[x].R = internal[y].R;			//
				internal[y].R = temp;				//
				internal[x].R->M = &internal[x];	// - adjust parent pointers
				internal[y].R->M = &internal[y];			//
				d->swapEdges(y, internal[y].L->index, LEFT, y,
						internal[y].R->index, RIGHT);
				temp = internal[y].L;			// - swap L and R of [y]
				internal[y].L = internal[y].R;			//
				internal[y].R = temp;				//
				internal[y].n = n_i + n_j;			// - update n
				internal[x].e = e_x;					// - update e_i
				internal[y].e = e_y;					//
				internal[x].p = p_x;					// - update p_i
				internal[y].p = p_y;					//
				internal[x].logL = L_x;					// - update logL_i
				internal[y].logL = L_y;					//
				internal[y].label = internal[x].label;		// - order-property
				L += dLogL;				// - update LogL
				delta = dLogL;				//

				// TRAP: Catches violations of the ordering property
				if (internal[x].label > internal[x].L->label
						|| internal[y].label > internal[y].L->label) {
					printDendrogram();
					if (internal[x].label > internal[x].L->label) {
						std::cout
								<< "**** WARNING - order property violated by internal["
								<< x << "]" << std::endl;
						std::cout << "x    (p = " << internal[x].p << "\te = "
								<< internal[x].e << "\tnL = ";
						std::cout << internal[x].L->n << "\tnR = "
								<< internal[x].R->n << "\tlabel = "
								<< internal[x].label << ")\tinternal[" << x
								<< "]\t(D)" << std::endl;
						if (internal[x].L->type == GRAPH) {
							std::cout << "x->L [" << internal[x].L->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->L (p = " << internal[x].L->p
									<< "\te = " << internal[x].L->e
									<< "\tnL = ";
							std::cout << internal[x].L->L->n << "\tnR = "
									<< internal[x].L->R->n << "\tlabel = "
									<< internal[x].L->label << ")\tinternal["
									<< internal[x].L->index << "]\t(D)"
									<< std::endl;
						}
						if (internal[x].R->type == GRAPH) {
							std::cout << "x->R [" << internal[x].R->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->R (p = " << internal[x].R->p
									<< "\te = " << internal[x].R->e
									<< "\tnL = ";
							std::cout << internal[x].R->L->n << "\tnR = "
									<< internal[x].R->R->n << "\tlabel = "
									<< internal[x].R->label << ")\tinternal["
									<< internal[x].R->index << "]\t(D)"
									<< std::endl;
						}
					}
					if (internal[y].label > internal[y].L->label) {
						std::cout
								<< "**** WARNING - order property violated by internal["
								<< y << "]" << std::endl;
						std::cout << "y    (p = " << internal[y].p << "\te = "
								<< internal[y].e << "\tnL = ";
						std::cout << internal[y].L->n << "\tnR = "
								<< internal[y].R->n << "\tlabel = "
								<< internal[y].label << ")\tinternal[" << y
								<< "]\t(D)" << std::endl;
						if (internal[y].L->type == GRAPH) {
							std::cout << "y->L [" << internal[y].L->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->L (p = " << internal[y].L->p
									<< "\te = " << internal[y].L->e
									<< "\tnL = ";
							std::cout << internal[y].L->L->n << "\tnR = "
									<< internal[y].L->R->n << "\tlabel = "
									<< internal[y].L->label << ")\tinternal["
									<< internal[y].L->index << "]\t(D)"
									<< std::endl;
						}
						if (internal[y].R->type == GRAPH) {
							std::cout << "y->R [" << internal[y].R->index
									<< "]\t(G)" << std::endl;
						} else {
							std::cout << "i->R (p = " << internal[y].R->p
									<< "\te = " << internal[y].R->e
									<< "\tnL = ";
							std::cout << internal[y].R->L->n << "\tnR = "
									<< internal[y].R->R->n << "\tlabel = "
									<< internal[y].R->label << ")\tinternal["
									<< internal[y].R->index << "]\t(D)"
									<< std::endl;
						}
					}
					return false;
				}
			}
		}
	}

	return true;
}

// ********************************************************************************************************

void dendro::printDendrogram() {
	std::cout << "\nLEAFS = " << n << std::endl << "# ";
	printSubTree(root);
	return;
}

// ********************************************************************************************************

void dendro::printSubTree(elementd *z) {
	if (z != nullptr) {
		if (z->type == GRAPH) {
			std::cout << "[" << z->label << "]" << std::endl; //"\t(" << z->L << " " << z->R << ") - " << z->M <<  std::endl;
			return;
		} else if (z->type == DENDRO) {
			std::cout << "(p = " << z->p << "\te = " << z->e << "\tnL = "
					<< z->L->n << "\tnR = " << z->R->n << "\tlabel = "
					<< z->label << ")\tinternal[" << z->index << "]"
					<< std::endl; // "\t(" << z->L << " " << z->R << ") - " << z->M <<  std::endl;
			std::cout << "L ";
			printSubTree(z->L);
			std::cout << std::endl;
			std::cout << "R ";
			printSubTree(z->R);
			std::cout << std::endl;
		} else {
			std::cout << "(p = " << z->p << "\te = " << z->e << "\tnL = "
					<< z->L->n << "\tnR = " << z->R->n << "\tlabel = "
					<< z->label << ")\tinternal[" << z->index << "] " << z->type
					<< std::endl;
			std::cout << "L ";
			printSubTree(z->L);
			std::cout << std::endl;
			std::cout << "R ";
			printSubTree(z->R);
			std::cout << std::endl;
		}
	}
	return;
}

// ********************************************************************************************************

void dendro::QsortMain(block* array, int left, int right) {
	if (right > left) {
		int pivot = left;
		int part = QsortPartition(array, left, right, pivot);
		QsortMain(array, left, part - 1);
		QsortMain(array, part + 1, right);
	}
	return;
}

int dendro::QsortPartition(block* array, int left, int right, int index) {
	block p_value, temp;
	p_value.x = array[index].x;
	p_value.y = array[index].y;

	// swap(array[p_value], array[right])
	temp.x = array[right].x;
	temp.y = array[right].y;
	array[right].x = array[index].x;
	array[right].y = array[index].y;
	array[index].x = temp.x;
	array[index].y = temp.y;

	int stored = left;
	for (int i = left; i < right; i++) {
		if (array[i].x <= p_value.x) {
			// swap(array[stored], array[i])
			temp.x = array[i].x;
			temp.y = array[i].y;
			array[i].x = array[stored].x;
			array[i].y = array[stored].y;
			array[stored].x = temp.x;
			array[stored].y = temp.y;
			stored++;
		}
	}
	// swap(array[right], array[stored])
	temp.x = array[stored].x;
	temp.y = array[stored].y;
	array[stored].x = array[right].x;
	array[stored].y = array[right].y;
	array[right].x = temp.x;
	array[right].y = temp.y;

	return stored;
}

// ********************************************************************************************************

void dendro::refreshLikelihood() {
	// recalculates the log-likelihood of the dendrogram structure
	L = 0.0;
	double dL;
	int nL_nR, ei;
	for (int i = 0; i < (n - 1); i++) {
		nL_nR = internal[i].L->n * internal[i].R->n;
		ei = internal[i].e;
		internal[i].p = (double) (ei) / (double) (nL_nR);
		if (ei == 0 || ei == nL_nR) {
			dL = 0.0;
		} else {
			dL = ei * log(internal[i].p)
					+ (nL_nR - ei) * log(1.0 - internal[i].p);
		}
		internal[i].logL = dL;
		L += dL;
	}
	return;
}

// ********************************************************************************************************

void dendro::resetDendrograph() {
	// Reset the dendrograph structure for the next trial
	if (leaf != nullptr) {
		delete[] leaf;
		leaf = nullptr;
	}    // O(n)
	if (internal != nullptr) {
		delete[] internal;
		internal = nullptr;
	}    // O(n)
	if (d != nullptr) {
		delete d;
		d = nullptr;
	}    // O(n)
	root = nullptr;
	if (paths != nullptr) {
		list *curr, *prev;
		for (int i = 0; i < n; i++) {
			curr = paths[i];
			while (curr != nullptr) {
				prev = curr;
				curr = curr->next;
				delete prev;
				prev = nullptr;
			}
			paths[i] = nullptr;
		}
		delete[] paths;
	}
	paths = nullptr;
	L = 1.0;

	return;
}

// ********************************************************************************************************

list* dendro::reversePathToRoot(const int leafIndex) {
	list *head, *subhead, *newlist;
	head = subhead = newlist = nullptr;
	elementd *current = &leaf[leafIndex];

	while (current != nullptr) {				// continue until we're finished
		newlist = new list;			// add this node to the path
		newlist->x = current->index;
		newlist->next = nullptr;
		if (head == nullptr) {
			head = newlist;
		} else {
			subhead = head;
			head = newlist;
			head->next = subhead;
		}
		current = current->M;
	}
	return head;
}

// ********************************************************************************************************

void dendro::sampleAdjacencyLikelihoods() {
	// Here, we sample the probability values associated with every adjacency in A, weighted by
	// their likelihood. The weighted histogram is stored in the graph data structure, so we simply
	// need to add an observation to each node-pair that corresponds to the associated branch point's
	// probability and the dendrogram's overall likelihood.
	bool flag_debug = false;
//	int num_bins = g->getNumBins();
	double nn;
	double norm = ((double) (n) * (double) (n)) / 4.0;
//	char pauseme;

	if (flag_debug) {
		std::cout << "dendro:: tabulating A'" << std::endl;
	}
	if (L > 0.0) {
		L = 0.0;
	}
	elementd* ancestor;
	list *currL, *prevL;
	if (paths != nullptr) {
		for (int i = 0; i < n; i++) {
			currL = paths[i];
			while (currL != nullptr) {
				prevL = currL;
				currL = currL->next;
				delete prevL;
				prevL = nullptr;
			}
			paths[i] = nullptr;
		}
		delete[] paths;
	}
	paths = nullptr;
	paths = new list*[n];
	for (int i = 0; i < n; i++) {
		paths[i] = reversePathToRoot(i);
	}	// construct paths from root, O(n^2) at worst
	for (int i = 0; i < n; i++) {// add obs for every node-pair, always O(n^2)
		for (int j = i + 1; j < n; j++) {
			ancestor = findCommonAncestor(paths, i, j);	// find internal node, O(n) at worst
			nn = ((double) (ancestor->L->n) * (double) (ancestor->R->n)) / norm;
			g->addAdjacencyObs(i, j, ancestor->p, nn);// add obs of ->p to (i,j) histogram, and
			g->addAdjacencyObs(j, i, ancestor->p, nn);// add obs of ->p to (j,i) histogram
		}												//
	}													//
	g->addAdjacencyEnd();		// finish-up: upate total weight in histograms

	return;
}

