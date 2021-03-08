/*
 recommend.c
 $LastChangedDate: 2008-10-22 17:39:18 -0500 (Wed, 22 Oct 2008) $
 $Revision: 134 $
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <search.h>
#include <stdexcept>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>

#include "linkpred/predictors/undirected/usbmpredictor/tools.h"
#include "linkpred/predictors/undirected/usbmpredictor/graph.h"
#include "linkpred/predictors/undirected/usbmpredictor/modules.h"
#include "linkpred/predictors/undirected/usbmpredictor/models.h"
#include "linkpred/predictors/undirected/usbmpredictor/missing.h"
#include "linkpred/predictors/undirected/usbmpredictor/recommend.h"

#define max(A, B) ((A > B)? A : B)

/*
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 Auxiliary functions
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 */
struct query *
CreateQuery(struct node_gra *node1, struct node_gra *node2) {
	struct query *theQuery;
	theQuery = (struct query *) calloc(1, sizeof(struct query));
	theQuery->n1 = node1;
	theQuery->n2 = node2;
	theQuery->score = -1.;
	return theQuery;
}

void FreeQuery(struct query *q) {
	free(q);
	return;
}

/*
 -----------------------------------------------------------------------------
 Given a (bipartite) network of ratings, return the number of
 unobserved pairs.
 -----------------------------------------------------------------------------
 */
int CountUnobserved(struct binet *ratings) {
	return (CountNodes(ratings->net1) * CountNodes(ratings->net2))
			- NLinksBipart(ratings);
}

/*
 -----------------------------------------------------------------------------
 Given a (bipartite) network of ratings, return a list of unobserved
 pairs.
 -----------------------------------------------------------------------------
 */
struct query **
BuildUnobservedSet(struct binet *ratings) {
	struct node_gra *p, *m;
	int nunobserved = CountUnobserved(ratings), n = 0;
	struct query **unobservedSet = NULL;

	unobservedSet = (struct query **) calloc(nunobserved,
			sizeof(struct query *));
	p = ratings->net1;
	while ((p = p->next) != NULL) {
		m = ratings->net2;
		while ((m = m->next) != NULL) {
			if (IsThereLink(p, m) == 0) {
				unobservedSet[n++] = CreateQuery(p, m);
			}
		}
	}
	return unobservedSet;
}

/*
 -----------------------------------------------------------------------------
 Given a (bipartite) network of ratings, remove all ratings whose
 value is r.
 -----------------------------------------------------------------------------
 */
void RemoveRatings(struct binet *ratings, int r) {
	struct node_gra *p = NULL;
	struct node_lis *m = NULL;

	p = ratings->net1;
	while ((p = p->next) != NULL) {
		m = p->neig;
		while (m->next != NULL) {
			if ((m->next)->weight == r) {
				RemoveLink(p, (m->next)->ref, 1);
			} else {
				m = m->next;
			}
		}
	}
	return;
}

/*
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 I/O functions
 ---------------------------------------------------------------------
 ---------------------------------------------------------------------
 */
/*
 -----------------------------------------------------------------------------
 From a file with format:

 p1 m1 0
 p1 m2 5
 ...

 (representing person 1 rating movie 1 with a 0, movie 2 with a 5,
 and so on) return a bipartite network with link weights representing
 the ratings.
 -----------------------------------------------------------------------------
 */
struct binet *
ReadRecommenderObservations(FILE *inFile) {
	return FBuildNetworkBipart(inFile, 1, 0);
}

/*
 -----------------------------------------------------------------------------
 Read a file with a query list with format:

 p1 m1
 p1 m2
 p2 m3
 ...

 representing the pairs (person, movie) that we want to predict. If a
 node in a query is missing from the bipartite ratings network, it is
 added with no links.
 
 -----------------------------------------------------------------------------
 */
struct query **
ReadQueries(FILE *inFile, int nQueries, struct binet *binet) {
	void *dict1, *dict2;
	char node1[MAX_LABEL_LENGTH], node2[MAX_LABEL_LENGTH];
	struct query **querySet;
//	struct query *q = NULL;
	int nq;
	struct node_gra *n1 = NULL, *n2 = NULL;

	/* Create the search dictionaries */
	dict1 = MakeLabelDict(binet->net1);
	dict2 = MakeLabelDict(binet->net2);

	/* Allocate space for the queries */
	querySet = (struct query **) calloc(nQueries, sizeof(struct query *));

	/* Go through the file and create the queries */
	for (nq = 0; nq < nQueries; nq++) {
		if (fscanf(inFile, "%s %s\n", &node1[0], &node2[0]) == EOF) {
			throw std::runtime_error("Could not read from file");
		}
		n1 = GetNodeDict(node1, dict1);
		/* add to binet if absent */
		if (n1 == NULL) {
			n1 = CreateNodeGraph(binet->net1, node1);
			fprintf(stderr,
					"WARNING! query node %s in net1 was never observed: adding!\n",
					n1->label);
			FreeLabelDict(dict1);
			dict1 = MakeLabelDict(binet->net1);
		}
		n2 = GetNodeDict(node2, dict2);
		/* add to binet if absent */
		if (n2 == NULL) {
			n2 = CreateNodeGraph(binet->net2, node2);
			fprintf(stderr,
					"WARNING! query node %s in net2 was never observed: adding!\n",
					n2->label);
			FreeLabelDict(dict2);
			dict2 = MakeLabelDict(binet->net2);
		}
		querySet[nq] = CreateQuery(n1, n2);
	}

	/* Free memory */
	FreeLabelDict(dict1);
	FreeLabelDict(dict2);

	/* Done */
	return querySet;
}

/*
 -----------------------------------------------------------------------------
 -----------------------------------------------------------------------------
 Recommender functions
 -----------------------------------------------------------------------------
 -----------------------------------------------------------------------------
 -----------------------------------------------------------------------------
 Hamiltoninan H for 2-state recommender. The part1 and part2
 partitions must be mapped onto a bipartite network whose links
 represent observed ratings. The weight of the link represents the
 rating. CAUTION!! ALL LINKS IN THE RATINGS BIPARTITE NETWORK ARE
 CONSIDERED OBSERVED--IF SOME ARE IN THE QUERY SET, THEY SHOULD BE
 REMOVED BEFORE GETTING HERE!!!
 -----------------------------------------------------------------------------
 */
double H2State(struct group *part1, struct group *part2) {
	struct group *g1 = NULL, *g2 = NULL;
	double r, l, H = 0.0;
	int ng1 = 0, ng2 = 0; /* Number of non-empty groups */
	int nnod1 = 0, nnod2 = 0; /* Number of nodes */

	/* Go through all group pairs */
	g1 = part1;
	while ((g1 = g1->next) != NULL) {
		if (g1->size > 0) {
			ng1++;
			nnod1 += g1->size;
			g2 = part2;
			while ((g2 = g2->next) != NULL) {
				if (g2->size > 0) {
					r = NG2GLinks(g1, g2); /* The total number of observations */
					l = NWeightG2GLinks(g1, g2, (double) 1); /* The 1s */
					H += log(r + 1) + LogChoose(r, l);
				}
			}
		}
	}

	/* Add the labeled-group sampling corrections */
	g2 = part2;
	while ((g2 = g2->next) != NULL) {
		if (g2->size > 0) {
			ng2++;
			nnod2 += g2->size;
		}
	}
	H -= gsl_sf_lnfact(nnod1 - ng1) + gsl_sf_lnfact(nnod2 - ng2);

	/* Done */
	return H;
}

/*
 -----------------------------------------------------------------------------
 Do a Monte Carlo step for the 2-state recommender system. Fast
 version for sparse observation matrices.
 -----------------------------------------------------------------------------
 */
void MCStep2State(int factor, double *H, struct node_gra **nlist1,
		struct node_gra **nlist2, struct group **glist1, struct group **glist2,
		struct group *part1, struct group *part2, int nnod1, int nnod2,
		int *ng1, int *ng2, int **N1G2_0, int **N2G1_0, int **N1G2_1,
		int **N2G1_1, int **G1G2_0, int **G2G1_0, int **G1G2_1, int **G2G1_1,
		double *LogList, int LogListSize, double **LogChooseList,
		int LogChooseListSize, double *LogFactList, int LogFactListSize,
		gsl_rng *gen) {
	double dH;
	struct group *oldg, *newg, *g, *g2;
	int ***N2G_0, ***N2Ginv_0, ***N2G_1, ***N2Ginv_1;
	int ***G2G_0, ***G2Ginv_0, ***G2G_1, ***G2Ginv_1;
	int move;
	double set_ratio;
	struct node_gra *node = NULL;
	struct node_lis *nei = NULL;
	int dice, r, l;
	int oldgnum, newgnum; // q;
	int i, nnod, ngroup;
	int *ng; /* Number of non-empty groups (for labeled-group correction) */
//	int j;
//	int move_in;

	/* Ratio of moves in each of the sets */
	set_ratio = (double) (nnod1 * nnod1 - 1)
			/ (double) (nnod1 * nnod1 + nnod2 * nnod2 - 2);

	/* The moves */
	for (move = 0; move < (nnod1 + nnod2) * factor; move++) {

		/* The move */
		if (gsl_rng_uniform(gen) < set_ratio) { /* move in first set */
//			move_in = 1;
			G2G_0 = &G1G2_0; /* Group-to-group links */
			G2Ginv_0 = &G2G1_0;
			G2G_1 = &G1G2_1;
			G2Ginv_1 = &G2G1_1;
			N2G_0 = &N1G2_0; /* Node-to-group links */
			N2Ginv_0 = &N2G1_0;
			N2G_1 = &N1G2_1;
			N2Ginv_1 = &N2G1_1;
			nnod = nnod1;
			ng = ng1;
			ngroup = nnod2;
			dice = floor(gsl_rng_uniform(gen) * (double) nnod1);
			node = nlist1[dice];
			oldgnum = node->inGroup;
			do {
				newgnum = floor(gsl_rng_uniform(gen) * (double) nnod1);
			} while (newgnum == oldgnum);
			oldg = glist1[oldgnum];
			newg = glist1[newgnum];
			g = g2 = part2;
		} else { /* move in second set */
//			move_in = 2;
			G2G_0 = &G2G1_0; /* Group-to-group links */
			G2Ginv_0 = &G1G2_0;
			G2G_1 = &G2G1_1;
			G2Ginv_1 = &G1G2_1;
			N2G_0 = &N2G1_0; /* Node-to-group links */
			N2Ginv_0 = &N1G2_0;
			N2G_1 = &N2G1_1;
			N2Ginv_1 = &N1G2_1;
			nnod = nnod2;
			ng = ng2;
			ngroup = nnod1;
			dice = floor(gsl_rng_uniform(gen) * (double) nnod2);
			node = nlist2[dice];
			oldgnum = node->inGroup;
			do {
				newgnum = floor(gsl_rng_uniform(gen) * (double) nnod2);
			} while (newgnum == oldgnum);
			oldg = glist2[oldgnum];
			newg = glist2[newgnum];
			g = g2 = part1;
		}

		/* THE CHANGE OF ENERGY: OLD CONFIGURATION CONTRIBUTION */
		dH = 0.0;

		while ((g = g->next) != NULL) {
			if (g->size > 0) { /* group is not empty */
				/* old configuration, old group */
				r = (*G2G_0)[oldgnum][g->label] + (*G2G_1)[oldgnum][g->label];
				l = (*G2G_1)[oldgnum][g->label];
				dH -= FastLog(r + 1, LogList, LogListSize)
						+ FastLogChoose(r, l, LogChooseList, LogChooseListSize);
				/* old configuration, new group */
				r = (*G2G_0)[newgnum][g->label] + (*G2G_1)[newgnum][g->label];
				l = (*G2G_1)[newgnum][g->label];
				dH -= FastLog(r + 1, LogList, LogListSize)
						+ FastLogChoose(r, l, LogChooseList, LogChooseListSize);
			}
		}
		/* labeled-groups sampling correction */
		if ((oldg->size == 1) || (newg->size == 0)) {
			dH += FastLogFact(nnod - *ng, LogFactList, LogFactListSize);
		}

		/* Tentatively move the node to the new group and update G2G and
		 N2Ginv matrices */
		MoveNode(node, oldg, newg);
		for (i = 0; i < ngroup; i++) { /* update G2G links */
			(*G2G_0)[oldgnum][i] -= (*N2G_0)[node->num][i];
			(*G2G_0)[newgnum][i] += (*N2G_0)[node->num][i];
			(*G2Ginv_0)[i][oldgnum] -= (*N2G_0)[node->num][i];
			(*G2Ginv_0)[i][newgnum] += (*N2G_0)[node->num][i];
			(*G2G_1)[oldgnum][i] -= (*N2G_1)[node->num][i];
			(*G2G_1)[newgnum][i] += (*N2G_1)[node->num][i];
			(*G2Ginv_1)[i][oldgnum] -= (*N2G_1)[node->num][i];
			(*G2Ginv_1)[i][newgnum] += (*N2G_1)[node->num][i];
		}
		nei = node->neig; /* update N2Ginv links */
		while ((nei = nei->next) != NULL) {
			if (nei->weight == (double) 0) {
				(*N2Ginv_0)[nei->ref->num][oldgnum] -= 1;
				(*N2Ginv_0)[nei->ref->num][newgnum] += 1;
			} else if (nei->weight == (double) 1) {
				(*N2Ginv_1)[nei->ref->num][oldgnum] -= 1;
				(*N2Ginv_1)[nei->ref->num][newgnum] += 1;
			} else {
				fprintf(stderr, "ERROR!!!");
			}
		}
		if (oldg->size == 0) /* update number of non-empty groups */
			(*ng) -= 1;
		if (newg->size == 1)
			(*ng) += 1;

		/* THE CHANGE OF ENERGY: NEW CONFIGURATION CONTRIBUTION */

		while ((g2 = g2->next) != NULL) {
			if (g2->size > 0) { /* group is not empty */
				/* new configuration, old group */
				r = (*G2G_0)[oldgnum][g2->label] + (*G2G_1)[oldgnum][g2->label];
				l = (*G2G_1)[oldgnum][g2->label];
				dH += FastLog(r + 1, LogList, LogListSize)
						+ FastLogChoose(r, l, LogChooseList, LogChooseListSize);
				/* new configuration, new group */
				r = (*G2G_0)[newgnum][g2->label] + (*G2G_1)[newgnum][g2->label];
				l = (*G2G_1)[newgnum][g2->label];
				dH += FastLog(r + 1, LogList, LogListSize)
						+ FastLogChoose(r, l, LogChooseList, LogChooseListSize);
			}
		}
		/* labeled-groups sampling correction */
		if ((oldg->size == 0) || (newg->size == 1)) {
			dH -= FastLogFact(nnod - *ng, LogFactList, LogFactListSize);
		}

		/* METROPOLIS ACCEPTANCE */
		if ((dH <= 0.0) || (gsl_rng_uniform(gen) < exp(-dH))) {
			/* accept move: update energy */
			*H += dH;
		} else {
			/* undo the move */
			MoveNode(node, newg, oldg);
			for (i = 0; i < ngroup; i++) { /* update G2G links */
				(*G2G_0)[oldgnum][i] += (*N2G_0)[node->num][i];
				(*G2G_0)[newgnum][i] -= (*N2G_0)[node->num][i];
				(*G2Ginv_0)[i][oldgnum] += (*N2G_0)[node->num][i];
				(*G2Ginv_0)[i][newgnum] -= (*N2G_0)[node->num][i];
				(*G2G_1)[oldgnum][i] += (*N2G_1)[node->num][i];
				(*G2G_1)[newgnum][i] -= (*N2G_1)[node->num][i];
				(*G2Ginv_1)[i][oldgnum] += (*N2G_1)[node->num][i];
				(*G2Ginv_1)[i][newgnum] -= (*N2G_1)[node->num][i];
			}
			nei = node->neig; /* update N2Ginv links */
			while ((nei = nei->next) != NULL) {
				if (nei->weight == (double) 0) {
					(*N2Ginv_0)[nei->ref->num][oldgnum] += 1;
					(*N2Ginv_0)[nei->ref->num][newgnum] -= 1;
				} else if (nei->weight == (double) 1) {
					(*N2Ginv_1)[nei->ref->num][oldgnum] += 1;
					(*N2Ginv_1)[nei->ref->num][newgnum] -= 1;
				}
			}
			if (oldg->size == 1) /* update number of non-empty groups */
				(*ng) += 1;
			if (newg->size == 0)
				(*ng) -= 1;
		}
	} /* Moves completed: done! */

	return;
}

/*
 ---------------------------------------------------------------------
 Get the decorrelation step necessary to sample decorrelated
 partitions. From part1 and part2, the longest decorrelation step is
 chosen.
 ---------------------------------------------------------------------
 */
int GetDecorrelationStep2State(double *H, struct node_gra **nlist1,
		struct node_gra **nlist2, struct group **glist1, struct group **glist2,
		struct group *part1, struct group *part2, int nnod1, int nnod2,
		int *ng1, int *ng2, int **N1G2_0, int **N2G1_0, int **N1G2_1,
		int **N2G1_1, int **G1G2_0, int **G2G1_0, int **G1G2_1, int **G2G1_1,
		double *LogList, int LogListSize, double **LogChooseList,
		int LogChooseListSize, double *LogFactList, int LogFactListSize,
		gsl_rng *gen, char verbose_sw) {
	struct group *part1Ref, *part2Ref;
	int step, x1, x2;
	double y11 = 0, y12 = 0, y21, y22;
//	double mutualInfo;
	int rep, nrep = 10;
	double *decay1, meanDecay1, *decay2, meanDecay2;
	double **decay, meanDecay, sigmaDecay, result;
	int norm = 0;

	x2 = (nnod1 + nnod2) / 5;
	if (x2 < 10)
		x2 = 10;
	x1 = x2 / 4;

	/* Get the nrep initial estimates */
	decay1 = allocate_d_vec(nrep);
	decay2 = allocate_d_vec(nrep);
	for (rep = 0; rep < nrep; rep++) {
		switch (verbose_sw) {
		case 'q':
			break;
		default:
			fprintf(stderr, "#\n# Estimating decorrelation time (%d/%d)\n",
					rep + 1, nrep);
			break;
		}
		part1Ref = CopyPartition(part1);
		part2Ref = CopyPartition(part2);
		for (step = 0; step <= x2; step++) {
			switch (verbose_sw) {
			case 'd':
				fprintf(stderr, "# %d / %d\n", step, x2);
				break;
			default:
				break;
			}
			MCStep2State(1, H, nlist1, nlist2, glist1, glist2, part1, part2,
					nnod1, nnod2, ng1, ng2, N1G2_0, N2G1_0, N1G2_1, N2G1_1,
					G1G2_0, G2G1_0, G1G2_1, G2G1_1, LogList, LogListSize,
					LogChooseList, LogChooseListSize, LogFactList,
					LogFactListSize, gen);
			if (step == x1) {
				y11 = MutualInformation(part1Ref, part1);
				y12 = MutualInformation(part2Ref, part2);
			}
		}
		y21 = MutualInformation(part1Ref, part1);
		y22 = MutualInformation(part2Ref, part2);
		if (nnod1 > 1)
			decay1[rep] = 2. * CalculateDecay(nnod1, x1, y11, x2, y21);
		else
			decay1[rep] = 1.e-6;
		if (nnod2 > 1)
			decay2[rep] = 2. * CalculateDecay(nnod2, x1, y12, x2, y22);
		else
			decay2[rep] = 1.e-6;
		switch (verbose_sw) {
		case 'q':
			break;
		default:
			fprintf(stderr, "# Decorrelation times (estimate %d) = %g %g\n",
					rep + 1, decay1[rep], decay2[rep]);
			break;
		}
		if (decay1[rep] < 0. || decay2[rep] < 0.) {
			rep--;
			switch (verbose_sw) {
			case 'q':
				break;
			default:
				fprintf(stderr, "#\tignoring...\n");
				break;
			}
		}
		/* Free memory */
		RemovePartition(part1Ref);
		RemovePartition(part2Ref);
	}

	/* Get rid of bad estimates (Chauvenet criterion)  */
	meanDecay1 = mean(decay1, nrep);
	meanDecay2 = mean(decay2, nrep);
	if (meanDecay1 > meanDecay2) {
		meanDecay = meanDecay1;
		sigmaDecay = stddev(decay1, nrep);
		decay = &decay1;
		switch (verbose_sw) {
		case 'q':
			break;
		default:
			fprintf(stderr, "# Considering decorrelation from partition 1\n");
			break;
		}
	} else {
		meanDecay = meanDecay2;
		sigmaDecay = stddev(decay2, nrep);
		decay = &decay2;
		switch (verbose_sw) {
		case 'q':
			break;
		default:
			fprintf(stderr, "# Considering decorrelation from partition 2\n");
			break;
		}
	}
	result = meanDecay * nrep;
	for (rep = 0; rep < nrep; rep++) {
		if (fabs((*decay)[rep] - meanDecay) / sigmaDecay > 2) {
			result -= (*decay)[rep];
			switch (verbose_sw) {
			case 'q':
				break;
			default:
				fprintf(stderr, "# Disregarding estimate %d\n", rep + 1);
				break;
			}
		} else {
			norm++;
		}
	}

	switch (verbose_sw) {
	case 'q':
		break;
	default:
		fprintf(stderr, "# Decorrelation step: %d\n",
				(int) (result / norm + 0.5));
		break;
	}

	/* Clean up */
	free_d_vec(decay1);
	free_d_vec(decay2);

	return (int) (result / norm + 0.5);
}

/*
 ---------------------------------------------------------------------
 
 ---------------------------------------------------------------------
 */
void ThermalizeMC2State(int decorStep, double *H, struct node_gra **nlist1,
		struct node_gra **nlist2, struct group **glist1, struct group **glist2,
		struct group *part1, struct group *part2, int nnod1, int nnod2,
		int *ng1, int *ng2, int **N1G2_0, int **N2G1_0, int **N1G2_1,
		int **N2G1_1, int **G1G2_0, int **G2G1_0, int **G1G2_1, int **G2G1_1,
		double *LogList, int LogListSize, double **LogChooseList,
		int LogChooseListSize, double *LogFactList, int LogFactListSize,
		gsl_rng *gen, char verbose_sw) {
	double HMean0 = 1.e10, HStd0 = 1.e-10, HMean1, HStd1, *Hvalues;
	int rep, nrep = 20;
	int equilibrated = 0;

	Hvalues = allocate_d_vec(nrep);

	do {

		/* MC steps */
		for (rep = 0; rep < nrep; rep++) {
			MCStep2State(decorStep, H, nlist1, nlist2, glist1, glist2, part1,
					part2, nnod1, nnod2, ng1, ng2, N1G2_0, N2G1_0, N1G2_1,
					N2G1_1, G1G2_0, G2G1_0, G1G2_1, G2G1_1, LogList,
					LogListSize, LogChooseList, LogChooseListSize, LogFactList,
					LogFactListSize, gen);
			switch (verbose_sw) {
			case 'q':
				break;
			default:
				fprintf(stderr, "%lf\n", *H);
				break;
			}
			Hvalues[rep] = *H;
		}

		/* Check for equilibration */
		HMean1 = mean(Hvalues, nrep);
		HStd1 = stddev(Hvalues, nrep);
		if (HMean0 - HStd0 / sqrt(nrep) < HMean1 + HStd1 / sqrt(nrep)) {
			equilibrated++;
			switch (verbose_sw) {
			case 'q':
				break;
			default:
				fprintf(stderr, "#\tequilibrated (%d/5) H=%lf\n", equilibrated,
						HMean1);
				break;
			}
		} else {
			switch (verbose_sw) {
			case 'q':
				break;
			default:
				fprintf(stderr, "#\tnot equilibrated yet H0=%g+-%g H1=%g+-%g\n",
						HMean0, HStd0 / sqrt(nrep), HMean1, HStd1 / sqrt(nrep));
				break;
			}
			HMean0 = HMean1;
			HStd0 = HStd1;
			equilibrated = 0;
		}

	} while (equilibrated < 5);

	/* Clean up */
	free_d_vec(Hvalues);

	return;
}

/*
 -----------------------------------------------------------------------------
 Return the score p(A_ij=1|A^O) of a collection {(i,j)} querySet of
 links, for a 2-state system. The ratings are a bipartite network
 with links (corresponding to observations) that have values 0 or 1.
 -----------------------------------------------------------------------------
 */
double *
MultiLinkScore2State(struct binet *ratings, struct query **querySet, int nquery,
		int nIter, gsl_rng *gen, char verbose_sw, int decorStep) {
	int nnod1 = CountNodes(ratings->net1), nnod2 = CountNodes(ratings->net2);
//	int nn1, nn2;
	struct node_gra *net1 = NULL, *net2 = NULL;
	struct group *part1 = NULL, *part2 = NULL;
	struct node_gra *p1 = NULL, *p2 = NULL; // *node = NULL;
	struct node_gra **nlist1 = NULL, **nlist2 = NULL;
	struct group **glist1 = NULL, **glist2 = NULL;
	struct group *lastg = NULL;
	double H;
	int iter;
	double *score;
	int i, j;
	int **N1G2_0 = NULL, **N2G1_0 = NULL;
	int **N1G2_1 = NULL, **N2G1_1 = NULL;
	int **G1G2_0 = NULL, **G2G1_0 = NULL;
	int **G1G2_1 = NULL, **G2G1_1 = NULL;
	int LogChooseListSize = 500;
	double **LogChooseList = InitializeFastLogChoose(LogChooseListSize);
	int LogListSize = 5000;
	double *LogList = InitializeFastLog(LogListSize);
	int LogFactListSize = 2000;
	double *LogFactList = InitializeFastLogFact(LogFactListSize);
//	struct node_lis *n1 = NULL, *n2 = NULL;
	double contrib;
	int norm = 0;
//	int dice;
	int r, l;
	void *dict1 = NULL, *dict2 = NULL;
	struct binet *ratingsClean = NULL;
	int q;
	int ng1, ng2;
	FILE *outfile = NULL;

	/*
	 PRELIMINARIES
	 */

	/* Create a ratings bipartite network that DOES NOT contain whatever
	 observations are in the query set. Map the queries to the new
	 network */
	ratingsClean = CopyBipart(ratings);
	dict1 = MakeLabelDict(ratingsClean->net1);
	dict2 = MakeLabelDict(ratingsClean->net2);
	for (q = 0; q < nquery; q++) {
		querySet[q]->n1 = GetNodeDict(querySet[q]->n1->label, dict1);
		querySet[q]->n2 = GetNodeDict(querySet[q]->n2->label, dict2);
		if (IsThereLink(querySet[q]->n1, querySet[q]->n2) == 1) {
			RemoveLink(querySet[q]->n1, querySet[q]->n2, 1);
		}
	}
	FreeLabelDict(dict1);
	FreeLabelDict(dict2);

	/* Initialize scores */
	score = allocate_d_vec(nquery);
	for (q = 0; q < nquery; q++)
		score[q] = 0.0;

	/* Map nodes and groups to a list for faster access */
	fprintf(stderr, ">> Mapping nodes and groups to lists...\n");
	nlist1 = (struct node_gra **) calloc(nnod1, sizeof(struct node_gra *));
	glist1 = (struct group **) calloc(nnod1, sizeof(struct group *));
	lastg = part1 = CreateHeaderGroup();
	p1 = net1 = ratingsClean->net1;
	while ((p1 = p1->next) != NULL) {
		nlist1[p1->num] = p1;
		lastg = glist1[p1->num] = CreateGroup(lastg, p1->num);
	}
	nlist2 = (struct node_gra **) calloc(nnod2, sizeof(struct node_gra *));
	glist2 = (struct group **) calloc(nnod2, sizeof(struct group *));
	lastg = part2 = CreateHeaderGroup();
	p2 = net2 = ratingsClean->net2;
	while ((p2 = p2->next) != NULL) {
		nlist2[p2->num] = p2;
		lastg = glist2[p2->num] = CreateGroup(lastg, p2->num);
	}

	/* Place nodes in random partitions */
	fprintf(stderr, ">> Placing nodes in random partitions...\n");
	p1 = net1;
	ResetNetGroup(net1);
	while ((p1 = p1->next) != NULL) {
		/*     dice = floor(gsl_rng_uniform(gen) * (double)nnod1); */
		/*     AddNodeToGroup(glist1[dice], p1); */
		AddNodeToGroup(glist1[p1->num], p1);

	}
	p2 = net2;
	ResetNetGroup(net2);
	while ((p2 = p2->next) != NULL) {
		/*     dice = floor(gsl_rng_uniform(gen) * (double)nnod2); */
		/*     AddNodeToGroup(glist2[dice], p2); */
		AddNodeToGroup(glist2[p2->num], p2);
	}

	/* Get the initial group-to-group links matrix */
	fprintf(stderr, ">> Getting the initial group-to-group links matrix...\n");
	G1G2_0 = allocate_i_mat(nnod1, nnod2);
	G2G1_0 = allocate_i_mat(nnod2, nnod1);
	G1G2_1 = allocate_i_mat(nnod1, nnod2);
	G2G1_1 = allocate_i_mat(nnod2, nnod1);
	for (i = 0; i < nnod1; i++) {
		for (j = 0; j < nnod2; j++) {
			G1G2_0[i][j] = G2G1_0[j][i] = NWeightG2GLinks(glist1[i], glist2[j],
					(double) 0);
			G1G2_1[i][j] = G2G1_1[j][i] = NWeightG2GLinks(glist1[i], glist2[j],
					(double) 1);
		}
	}

	/* Get the initial node-to-group links matrix */
	fprintf(stderr, ">> Getting the initial node-to-group links matrix...\n");
	N1G2_0 = allocate_i_mat(nnod1, nnod2);
	N2G1_0 = allocate_i_mat(nnod2, nnod1);
	N1G2_1 = allocate_i_mat(nnod1, nnod2);
	N2G1_1 = allocate_i_mat(nnod2, nnod1);
	for (i = 0; i < nnod1; i++) {
		for (j = 0; j < nnod2; j++) {
			N1G2_0[i][j] = NWeightLinksToGroup(nlist1[i], glist2[j],
					(double) 0);
			N2G1_0[j][i] = NWeightLinksToGroup(nlist2[j], glist1[i],
					(double) 0);
			N1G2_1[i][j] = NWeightLinksToGroup(nlist1[i], glist2[j],
					(double) 1);
			N2G1_1[j][i] = NWeightLinksToGroup(nlist2[j], glist1[i],
					(double) 1);
		}
	}

	/* Get the initial number of non-empty groups */
	ng1 = NNonEmptyGroups(part1);
	ng2 = NNonEmptyGroups(part2);

	/*
	 GET READY FOR THE SAMPLING
	 */
	H = H2State(part1, part2);

	/* Get the decorrelation time */
	switch (verbose_sw) {
	case 'q':
		break;
	default:
		fprintf(stderr, "# CALCULATING DECORRELATION TIME\n");
		fprintf(stderr, "# ------------------------------\n");
		break;
	}
	if (decorStep <= 0) {
		decorStep = GetDecorrelationStep2State(&H, nlist1, nlist2, glist1,
				glist2, part1, part2, nnod1, nnod2, &ng1, &ng2, N1G2_0, N2G1_0,
				N1G2_1, N2G1_1, G1G2_0, G2G1_0, G1G2_1, G2G1_1, LogList,
				LogListSize, LogChooseList, LogChooseListSize, LogFactList,
				LogFactListSize, gen, verbose_sw);
	}

	/* Thermalization */
	switch (verbose_sw) {
	case 'q':
		break;
	default:
		fprintf(stderr, "#\n#\n# THERMALIZING\n");
		fprintf(stderr, "# ------------\n");
		break;
	}
	ThermalizeMC2State(decorStep, &H, nlist1, nlist2, glist1, glist2, part1,
			part2, nnod1, nnod2, &ng1, &ng2, N1G2_0, N2G1_0, N1G2_1, N2G1_1,
			G1G2_0, G2G1_0, G1G2_1, G2G1_1, LogList, LogListSize, LogChooseList,
			LogChooseListSize, LogFactList, LogFactListSize, gen, verbose_sw);

	/*
	 SAMPLIN' ALONG
	 */
	switch (verbose_sw) {
	case 'd':
		break;
	default:
		H = 0; /* Reset the origin of energies to avoid huge exponentials */
		break;
	}
	for (iter = 0; iter < nIter; iter++) {
		MCStep2State(decorStep, &H, nlist1, nlist2, glist1, glist2, part1,
				part2, nnod1, nnod2, &ng1, &ng2, N1G2_0, N2G1_0, N1G2_1, N2G1_1,
				G1G2_0, G2G1_0, G1G2_1, G2G1_1, LogList, LogListSize,
				LogChooseList, LogChooseListSize, LogFactList, LogFactListSize,
				gen);
		switch (verbose_sw) {
		case 'q':
			break;
		case 'v':
			fprintf(stderr, "%d %lf\n", iter, H);
			break;
		case 'd':
			fprintf(stderr, "%d %lf %lf\n", iter, H, H2State(part1, part2));
			FPrintPartition(stderr, part1, 0);
			FPrintPartition(stderr, part2, 0);
			break;
		}

		/* Check if the energy has gone below a certain threshold and, if
		 so, reset energies and start over */
		if (H < -400.0) {
			switch (verbose_sw) {
			case 'q':
				break;
			default:
				fprintf(stderr,
						"# System was not properly thermalized: starting over :(\n\n");
				break;
			}
			iter = 0;
			H = 0.0;
			norm = 0;
			for (q = 0; q < nquery; q++) {
				score[q] = 0.0;
			}
		}

		/* Update normalization */
		norm += 1;

		/* Update the scores */
		for (q = 0; q < nquery; q++) {
			l = G1G2_1[querySet[q]->n1->inGroup][querySet[q]->n2->inGroup];
			r =
					G1G2_0[querySet[q]->n1->inGroup][querySet[q]->n2->inGroup]
							+ G1G2_1[querySet[q]->n1->inGroup][querySet[q]->n2->inGroup];
			contrib = (float) (l + 1) / (float) (r + 2);
			score[q] += contrib;
		}

		/* Output temporary scores and partitions */
		if (iter % 100 == 0) {
			switch (verbose_sw) {
			case 'q':
				break;
			default:
				outfile = fopen("scores.tmp", "w");
				for (q = 0; q < nquery; q++) {
					fprintf(outfile, "%s %s %lf\n", ((querySet[q])->n1)->label,
							((querySet[q])->n2)->label,
							score[q] / (double) norm);
				}
				fclose(outfile);
				outfile = fopen("part1.tmp", "w");
				FPrintPartition(outfile, part1, 1);
				fclose(outfile);
				outfile = fopen("part2.tmp", "w");
				FPrintPartition(outfile, part2, 1);
				fclose(outfile);
				break;
			}
		}
	} /* End of iter loop */

	/* Normalize the scores */
	for (q = 0; q < nquery; q++)
		score[q] /= (double) norm;

	/* Remap the queries to the original network */
	dict1 = MakeLabelDict(ratings->net1);
	dict2 = MakeLabelDict(ratings->net2);
	for (q = 0; q < nquery; q++) {
		querySet[q]->n1 = GetNodeDict(querySet[q]->n1->label, dict1);
		querySet[q]->n2 = GetNodeDict(querySet[q]->n2->label, dict2);
	}
	FreeLabelDict(dict1);
	FreeLabelDict(dict2);

	/* Free dynamically allocated memory */
	RemovePartition(part1);
	RemovePartition(part2);
	free(glist1);
	free(glist2);
	free(nlist1);
	free(nlist2);
	free_i_mat(G1G2_0, nnod1);
	free_i_mat(G2G1_0, nnod2);
	free_i_mat(G1G2_1, nnod1);
	free_i_mat(G2G1_1, nnod2);
	free_i_mat(N1G2_0, nnod1);
	free_i_mat(N2G1_0, nnod2);
	free_i_mat(N1G2_1, nnod1);
	free_i_mat(N2G1_1, nnod2);
	FreeFastLog(LogList);
	FreeFastLogChoose(LogChooseList, LogChooseListSize);
	FreeFastLogFact(LogFactList);
	RemoveBipart(ratingsClean);

	/* Done */
	return score;
}
