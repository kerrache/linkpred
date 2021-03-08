/* sampling.h
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

#ifndef __SAMPLING_H__
#define __SAMPLING_H__

#include "linkpred/numerical/plfit/mt.hpp"
#include <stdlib.h>

namespace LinkPred {

namespace PLFit {
/**
 * Draws a sample from a uniform distribution over the [0; 1) interval.
 *
 * The interval is closed from the left and open from the right.
 *
 * \param  rng  the Mersenne Twister random number generator to use
 * \return the value drawn from the given uniform distribution.
 */
inline double plfit_runif_01(mt_rng_t* rng) {
	if (rng == 0) {
		return rand() / ((double) RAND_MAX);
	}
	return mt_uniform_01(rng);
}

/**
 * Draws a sample from a binomial distribution with the given count and
 * probability values.
 *
 * This function is borrowed from R; see the corresponding license in
 * \c rbinom.c. The return value is always an integer.
 *
 * The function is \em not thread-safe.
 *
 * \param  n    the number of trials
 * \param  p    the success probability of each trial
 * \param  rng  the Mersenne Twister random number generator to use
 * \return the value drawn from the given binomial distribution.
 */
double plfit_rbinom(double n, double p, mt_rng_t* rng);

/**
 * Draws a sample from a Pareto distribution with the given minimum value and
 * power-law exponent.
 *
 * \param  xmin    the minimum value of the distribution. Must be positive.
 * \param  alpha   the exponent. Must be positive
 * \param  rng     the Mersenne Twister random number generator to use
 *
 * \return the sample or NaN if one of the parameters is invalid
 */
inline double plfit_rpareto(double xmin, double alpha, mt_rng_t* rng) {
	if (alpha <= 0 || xmin <= 0)
		return NAN;

	/* 1-u is used in the base here because we want to avoid the case of
	 * sampling zero */
	return pow(1 - plfit_runif_01(rng), -1.0 / alpha) * xmin;
}

/**
 * Draws a given number of samples from a Pareto distribution with the given
 * minimum value and power-law exponent.
 *
 * \param  xmin    the minimum value of the distribution. Must be positive.
 * \param  alpha   the exponent. Must be positive
 * \param  n       the number of samples to draw
 * \param  rng     the Mersenne Twister random number generator to use
 * \param  result  the array where the result should be written. It must
 *                 have enough space to store n items
 *
 * \return \c PLFIT_EINVAL if one of the parameters is invalid, zero otherwise
 */
int plfit_rpareto_array(double xmin, double alpha, size_t n, mt_rng_t* rng,
		double* result);

/**
 * Draws a sample from a zeta distribution with the given minimum value and
 * power-law exponent.
 *
 * \param  xmin    the minimum value of the distribution. Must be positive.
 * \param  alpha   the exponent. Must be positive
 * \param  rng     the Mersenne Twister random number generator to use
 *
 * \return the sample or NaN if one of the parameters is invalid
 */
inline double plfit_rzeta(long int xmin, double alpha, mt_rng_t* rng) {
	double u, v, t;
	long int x;
	double alpha_minus_1 = alpha - 1;
	double minus_1_over_alpha_minus_1 = -1.0 / (alpha - 1);
	double b;
	double one_over_b_minus_1;

	if (alpha <= 0 || xmin < 1)
		return NAN;

	xmin = (long int) round(xmin);

	/* Rejection sampling for the win. We use Y=floor(U^{-1/alpha} * xmin) as the
	 * envelope distribution, similarly to Chapter X.6 of Luc Devroye's book
	 * (where xmin is assumed to be 1): http://luc.devroye.org/chapter_ten.pdf
	 *
	 * Some notes that should help me recover what I was doing:
	 *
	 * p_i = 1/zeta(alpha, xmin) * i^-alpha
	 * q_i = (xmin/i)^{alpha-1} - (xmin/(i+1))^{alpha-1}
	 *     = (i/xmin)^{1-alpha} - ((i+1)/xmin)^{1-alpha}
	 *     = [i^{1-alpha} - (i+1)^{1-alpha}] / xmin^{1-alpha}
	 *
	 * p_i / q_i attains its maximum at xmin=i, so the rejection constant is:
	 *
	 * c = p_xmin / q_xmin
	 *
	 * We have to accept the sample if V <= (p_i / q_i) * (q_xmin / p_xmin) =
	 * (i/xmin)^-alpha * [xmin^{1-alpha} - (xmin+1)^{1-alpha}] / [i^{1-alpha} - (i+1)^{1-alpha}] =
	 * [xmin - xmin^alpha / (xmin+1)^{alpha-1}] / [i - i^alpha / (i+1)^{alpha-1}] =
	 * xmin/i * [1-(xmin/(xmin+1))^{alpha-1}]/[1-(i/(i+1))^{alpha-1}]
	 *
	 * In other words (and substituting i with X, which is the same),
	 *
	 * V * (X/xmin) <= [1 - (1+1/xmin)^{1-alpha}] / [1 - (1+1/i)^{1-alpha}]
	 *
	 * Let b := (1+1/xmin)^{alpha-1} and let T := (1+1/i)^{alpha-1}. Then:
	 *
	 * V * (X/xmin) <= [(b-1)/b] / [(T-1)/T]
	 * V * (X/xmin) * (T-1) / (b-1) <= T / b
	 *
	 * which is the same as in Devroye's book, except for the X/xmin term, and
	 * the definition of b.
	 */
	b = pow(1 + 1.0 / xmin, alpha_minus_1);
	one_over_b_minus_1 = 1.0 / (b - 1);
	do {
		do {
			u = plfit_runif_01(rng);
			v = plfit_runif_01(rng);
			/* 1-u is used in the base here because we want to avoid the case of
			 * having zero in x */
			x = (long int) floor(pow(1 - u, minus_1_over_alpha_minus_1) * xmin);
		} while (x < xmin);
		t = pow((x + 1.0) / x, alpha_minus_1);
	} while (v * x * (t - 1) * one_over_b_minus_1 * b > t * xmin);

	return x;
}

/**
 * Draws a given number of samples from a zeta distribution with the given
 * minimum value and power-law exponent.
 *
 * \param  xmin    the minimum value of the distribution. Must be positive.
 * \param  alpha   the exponent. Must be positive
 * \param  n       the number of samples to draw
 * \param  rng     the Mersenne Twister random number generator to use
 * \param  result  the array where the result should be written. It must
 *                 have enough space to store n items
 *
 * \return \c PLFIT_EINVAL if one of the parameters is invalid, zero otherwise
 */
int plfit_rzeta_array(long int xmin, double alpha, size_t n, mt_rng_t* rng,
		double* result);

/**
 * Draws a sample from a uniform distribution with the given lower and
 * upper bounds.
 *
 * The lower bound is inclusive, the uppoer bound is not.
 *
 * \param  lo   the lower bound
 * \param  hi   the upper bound
 * \param  rng  the Mersenne Twister random number generator to use
 * \return the value drawn from the given uniform distribution.
 */
inline double plfit_runif(double lo, double hi, mt_rng_t* rng) {
	if (rng == 0) {
		return lo + rand() / ((double) RAND_MAX) * (hi - lo);
	}
	return lo + mt_uniform_01(rng) * (hi - lo);
}

/**
 * Random sampler using Walker's alias method.
 */
typedef struct {
	long int num_bins; /**< Number of bins */
	long int* indexes; /**< Index of the "other" element in each bin */
	double* probs; /**< Probability of drawing the "own" element from a bin */
} plfit_walker_alias_sampler_t;

/**
 * \brief Initializes the sampler with item probabilities.
 *
 * \param  sampler  the sampler to initialize
 * \param  ps   pointer to an array containing a value proportional to the
 *              sampling probability of each item in the set being sampled.
 * \param  n    the number of items in the array
 * \return error code
 */
int plfit_walker_alias_sampler_init(plfit_walker_alias_sampler_t* sampler,
		double* ps, size_t n);

/**
 * \brief Destroys an initialized sampler and frees the allocated memory.
 *
 * \param  sampler  the sampler to destroy
 */
void plfit_walker_alias_sampler_destroy(plfit_walker_alias_sampler_t* sampler);

/**
 * \brief Draws a given number of samples from the sampler and writes them
 *        to a given array.
 *
 * \param  sampler  the sampler to use
 * \param  xs       pointer to an array where the sampled items should be
 *                  written
 * \param  n        the number of samples to draw
 * \param  rng      the Mersenne Twister random number generator to use
 * \return error code
 */
int plfit_walker_alias_sampler_sample(
		const plfit_walker_alias_sampler_t* sampler, long int* xs, size_t n,
		mt_rng_t* rng);
} /* namespace PLFit */

} /* namespace LinkPred */

#endif

