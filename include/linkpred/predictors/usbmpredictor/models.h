/*
 models.h
 $LastChangedDate: 2010-04-28 16:05:59 +0200 (Wed, 28 Apr 2010) $
 $Revision: 219 $
 */

#ifndef RGRAPH_MODELS_H
#define RGRAPH_MODELS_H 1

#include <gsl/gsl_rng.h>

struct node_gra *EmptyGraph(int S);
struct node_gra *ERGraph(int S, double p, gsl_rng *gen);
struct node_gra *PAGraph(int S, int m, gsl_rng *gen);
struct node_gra *UndirectedBlockGraph(int ngroup, int *gsize, double **q,
		char output_sw, gsl_rng *gen);
struct node_gra *GirvanNewmanGraph(int ngroup, int gsize, double kin,
		double kout, char output_sw, gsl_rng *gen);

#endif /* !RGRAPH_MODELS_H */
