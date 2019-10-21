/* sampling.c
 *
 * Copyright (C) 2012 Tamas Nepusz
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <math.h>
#include "linkpred/numerical/plfit/error.h"
#include "linkpred/numerical/plfit/sampling.h"
#include "linkpred/numerical/plfit/platform.h"

namespace PLFit {
int plfit_rpareto_array(double xmin, double alpha, size_t n, mt_rng_t* rng,
		double* result) {
	double gamma;

	if (alpha <= 0 || xmin <= 0)
		return PLFIT_EINVAL;

	if (result == 0 || n == 0)
		return PLFIT_SUCCESS;

	gamma = -1.0 / alpha;
	while (n > 0) {
		/* 1-u is used in the base here because we want to avoid the case of
		 * sampling zero */
		*result = pow(1 - plfit_runif_01(rng), gamma) * xmin;
		result++;
		n--;
	}

	return PLFIT_SUCCESS;
}

int plfit_rzeta_array(long int xmin, double alpha, size_t n, mt_rng_t* rng,
		double* result) {
	double u, v, t;
	long int x;
	double alpha_minus_1 = alpha - 1;
	double minus_1_over_alpha_minus_1 = -1.0 / (alpha - 1);
	double b, one_over_b_minus_1;

	if (alpha <= 0 || xmin < 1)
		return PLFIT_EINVAL;

	if (result == 0 || n == 0)
		return PLFIT_SUCCESS;

	/* See the comments in plfit_rzeta for an explanation of the algorithm
	 * below. */
	xmin = (long int) round(xmin);
	b = pow(1 + 1.0 / xmin, alpha_minus_1);
	one_over_b_minus_1 = 1.0 / (b - 1);

	while (n > 0) {
		do {
			do {
				u = plfit_runif_01(rng);
				v = plfit_runif_01(rng);
				/* 1-u is used in the base here because we want to avoid the case of
				 * having zero in x */
				x = (long int) floor(
						pow(1 - u, minus_1_over_alpha_minus_1) * xmin);
			} while (x < xmin); /* handles overflow as well */
			t = pow((x + 1.0) / x, alpha_minus_1);
		} while (v * x * (t - 1) * one_over_b_minus_1 * b > t * xmin);
		*result = x;
		if (x < 0)
			abort();
		result++;
		n--;
	}

	return PLFIT_SUCCESS;
}

int plfit_walker_alias_sampler_init(plfit_walker_alias_sampler_t* sampler,
		double* ps, size_t n) {
	double *p, *p2, *ps_end;
	double sum;
	long int *short_sticks, *long_sticks;
	long int num_short_sticks, num_long_sticks;
	size_t i;

	sampler->num_bins = n;

	ps_end = ps + n;

	/* Initialize indexes and probs */
	sampler->indexes = (long int*) calloc(n, sizeof(long int));
	if (sampler->indexes == 0) {
		return PLFIT_ENOMEM;
	}
	sampler->probs = (double*) calloc(n, sizeof(double));
	if (sampler->probs == 0) {
		free(sampler->indexes);
		return PLFIT_ENOMEM;
	}

	/* Normalize the probability vector; count how many short and long sticks
	 * are there initially */
	for (sum = 0.0, p = ps; p != ps_end; p++) {
		sum += *p;
	}
	sum = n / sum;

	num_short_sticks = num_long_sticks = 0;
	for (p = ps, p2 = sampler->probs; p != ps_end; p++, p2++) {
		*p2 = *p * sum;
		if (*p2 < 1) {
			num_short_sticks++;
		} else if (*p2 > 1) {
			num_long_sticks++;
		}
	}

	/* Allocate space for short & long stick indexes */
	long_sticks = (long int*) calloc(num_long_sticks, sizeof(long int));
	if (long_sticks == 0) {
		free(sampler->probs);
		free(sampler->indexes);
		return PLFIT_ENOMEM;
	}
	short_sticks = (long int*) calloc(num_long_sticks, sizeof(long int));
	if (short_sticks == 0) {
		free(sampler->probs);
		free(sampler->indexes);
		free(long_sticks);
		return PLFIT_ENOMEM;
	}

	/* Initialize short_sticks and long_sticks */
	num_short_sticks = num_long_sticks = 0;
	for (i = 0, p = sampler->probs; i < n; i++, p++) {
		if (*p < 1) {
			short_sticks[num_short_sticks++] = i;
		} else if (*p > 1) {
			long_sticks[num_long_sticks++] = i;
		}
	}

	/* Prepare the index table */
	while (num_short_sticks && num_long_sticks) {
		long int short_index, long_index;
		short_index = short_sticks[--num_short_sticks];
		long_index = long_sticks[num_long_sticks - 1];
		sampler->indexes[short_index] = long_index;
		sampler->probs[long_index] = /* numerical stability */
		(sampler->probs[long_index] + sampler->probs[short_index]) - 1;
		if (sampler->probs[long_index] < 1) {
			short_sticks[num_short_sticks++] = long_index;
			num_long_sticks--;
		}
	}

	/* Fix numerical stability issues */
	while (num_long_sticks) {
		i = long_sticks[--num_long_sticks];
		sampler->probs[i] = 1;
	}
	while (num_short_sticks) {
		i = short_sticks[--num_short_sticks];
		sampler->probs[i] = 1;
	}

	return PLFIT_SUCCESS;
}

void plfit_walker_alias_sampler_destroy(plfit_walker_alias_sampler_t* sampler) {
	if (sampler->indexes) {
		free(sampler->indexes);
		sampler->indexes = 0;
	}
	if (sampler->probs) {
		free(sampler->probs);
		sampler->probs = 0;
	}
}

int plfit_walker_alias_sampler_sample(
		const plfit_walker_alias_sampler_t* sampler, long int *xs, size_t n,
		mt_rng_t* rng) {
	double u;
	long int j;
	long int *x;

	x = xs;

	if (rng == 0) {
		/* Using built-in RNG */
		while (n > 0) {
			u = rand() / ((double) RAND_MAX);
			j = rand() % sampler->num_bins;
			*x = (u < sampler->probs[j]) ? j : sampler->indexes[j];
			n--;
			x++;
		}
	} else {
		/* Using Mersenne Twister */
		while (n > 0) {
			u = mt_uniform_01(rng);
			j = mt_random(rng) % sampler->num_bins;
			*x = (u < sampler->probs[j]) ? j : sampler->indexes[j];
			n--;
			x++;
		}
	}

	return PLFIT_SUCCESS;
}

} /* namespace PLFit */
