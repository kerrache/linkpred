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

#include <linkpred/predictors/undirected/uhyppredictor.hpp>
#include "linkpred/utils/randomgen.hpp"
#include "linkpred/utils/log.hpp"
#ifdef LINKPRED_WITH_OPENMP
#include <omp.h>
#endif
#include <linkpred/utils/miscutils.hpp>
#include <fstream>
#include <cmath>
#include <algorithm>

namespace LinkPred {


template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UHYPPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::init() {
	double avgDeg;
	std::size_t minDeg;
	G = std::make_shared<Graph>();

	// Adding nodes
	for (NodeID i = 0; i < net->getNbNodes(); i++) {
		auto u = std::make_shared<Node>(i);
		G->addNewNode(u);
	}

	// Adding edges
	for (auto it = net->edgesBegin(); it != net->edgesEnd(); ++it) {
		auto vertex1 = Network::start(*it);
		auto vertex2 = Network::end(*it);

		auto u = G->findNode(vertex1);
		auto v = G->findNode(vertex2);

		u->addAdjNode(v);
		v->addAdjNode(u);
	}

	net->getDegStat(minDeg, maxDeg, avgDeg);
	m = minDeg;
	L = avgDeg / 2 - m;
	std::vector<double> degs;
	degs.resize(net->getNbNodes());
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (NodeID i = 0; i < net->getNbNodes(); i++) {
		degs[i] = std::max((std::size_t) 1, net->getDeg(i));
	}

	// Power law fitting
	gamma = Utils::plFit(degs).first;
	logger(logDebug, "gamma: " << gamma)

	if (seed == 0) {
		RandomGen rng;
		seed = rng.getInt();
	}
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UHYPPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::learn() {
	double N, t = 0, m_t = 0, R_t = 0;
	N = net->getNbNodes();
	std::shared_ptr < Node > u, v;

	double bar_k = 2 * (m + L);
	double Beta = 1 / (gamma - 1);

	auto nodes = G->getNodeList();

	auto model = std::make_shared<Graph>();

	//Embed all nodes with degree >= 1. You can change the 1 in the line below to any other value.
	//For example, changing it to 3, will embed only the nodes will degrees k >= 3.
	for (std::size_t i = maxDeg; i >= 1; i--) {
		for (std::size_t j = 0; j < nodes.size(); j++) {
			if (nodes[j]->getAdjNodeList().size() == i) {

				t++;
				u = nodes[j];
				double r_t = (2 / zeta) * log(t);
				u->setRadius(r_t);
				u->setInitRadius(r_t);
				model->addNewNode(u);
				auto nodes2 = model->getNodeList();

				for (std::size_t k = 0; k < nodes2.size(); k++) {
					v = nodes2[k];
					if (u != v) {
						v->setRadius(
								Beta * v->getInitRadius() + (1 - Beta) * r_t);
					}
				}

				double I_t = (1 / (1 - Beta)) * (1 - std::pow(t, -(1 - Beta)));
				double L_t = (2 * L * (1 - Beta)
						/ ((2 * Beta - 1)
								* std::pow((1 - std::pow(N, -(1 - Beta))), 2)))
						* (std::pow(N / t, (2 * Beta - 1)) - 1)
						* (1 - std::pow(t, -(1 - Beta)));

				m_t = m + L_t;

				if (t == 1) {
					u->setAngle(rng.getDouble(0, 2 * MathPI));
					logger(logDebug,
							std::to_string(t) + " " + std::to_string(u->getId())
									+ " " + std::to_string(u->getAngle()) + " "
									+ std::to_string(u->getRadius()))
					u->setFlag(1);
					continue;
				}

				//To have a more widespread node angular distribution we can ignore the links between the first few nodes for which m_t >= t-1.
				//This is OK since these nodes can have almost any angle with high probability. If you do not want this, comment the three lines of code below.
				if (m_t >= t - 1) {
					u->setFlag(1);
				}

				R_t = r_t
						- (2 / zeta)
								* log((2 * T / sin(T * MathPI)) * (I_t / m_t));

				if (R_t < 0) {
					logger(logInfo, "Negative R_t? Setting it to r_t\n")
					R_t = r_t;
				}

				u->set_R(R_t);

				double step = (double) 1 / (t), maxlog_L = -1000000000;
				double theta = 0;

				while (theta <= 2 * MathPI) {

					double theta_t = theta, log_L = 0;
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for private(v) reduction(+:log_L) if (parallel)
#endif
					for (std::size_t k = 0; k < nodes2.size(); k++) {

						v = nodes2[k];

						if (v != u) {

							double theta_s = v->getAngle();
							double r_s = v->getRadius();

							double dtheta = fabs(theta_t - theta_s);
							if (dtheta > MathPI) {
								dtheta = 2 * MathPI - dtheta;
							}

							double x_st;
							if (dtheta == 0) {
								x_st = fabs(r_t - r_s);
							} else {
								x_st =
										(1 / zeta)
												* acosh(
														(cosh(zeta * r_s)
																* cosh(
																		zeta
																				* r_t))
																- (sinh(
																		zeta
																				* r_s)
																		* sinh(
																				zeta
																						* r_t)
																		* cos(
																				dtheta)));
							}

							double P_st =
									1
											/ (1
													+ exp(
															(zeta / (2 * T))
																	* (x_st
																			- R_t)));

							//If flags are set on both nodes we ignore the links between them (see line 180).
							if ((u->getFlag() == 1) && (v->getFlag() == 1)) {
								log_L = log_L + log(1 - P_st);
								continue;
							}

							if (u->checkLink(v)) {
								log_L = log_L + log(P_st);
							} else {
								log_L = log_L + log(1 - P_st);
							}
						}
					}

					if (log_L >= maxlog_L) {
						u->setAngle(theta_t);
						maxlog_L = log_L;
					}

					theta = theta + step;
				}

				logger(logDebug,
						std::to_string(t) + " " + std::to_string(u->getId())
								+ " " + std::to_string(u->getAngle()) + " "
								+ std::to_string(u->getRadius()))
			}
			//Endif.

		}	        //Endfor (for all nodes with degree i).

		//The section of the code below runs correction steps as described in the paper.
		//Put this section into comments if you don't want to run correction steps.
		if ((i == 60) || (i == 40) || (i == 20) || (i == 10)) {

			logger(logDebug,
					"Running correction step for all nodes with degree greater than "
							+ std::to_string(i))

			//We repeat each correction step a number of times, equal here to the average node degree (fewer times can also be beneficial.)
			for (int round = 1; round <= bar_k; round++) {

				//cout<<"round "<<round<<"\n";

				auto nodes2 = model->getNodeList();

				for (std::size_t k = 0; k < nodes2.size(); k++) {

					u = nodes2[k];

					double step = (double) 1 / (t), maxlog_L = -1000000000;
					double theta = 0;

					while (theta <= 2 * MathPI) {

						double theta_u = theta, log_L = 0;

#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for private(v) reduction(+:log_L) if (parallel)
#endif
						for (std::size_t l = 0; l < nodes2.size(); l++) {

							v = nodes2[l];

							if (v != u) {

								double theta_v = v->getAngle();
								double dtheta = fabs(theta_v - theta_u);
								if (dtheta > MathPI) {
									dtheta = 2 * MathPI - dtheta;
								}

								double x_vu;
								double R_init;

								if (u->getInitRadius() > v->getInitRadius()) {//u came after v.

									double r_u = u->getInitRadius();
									double r_v = Beta * v->getInitRadius()
											+ (1 - Beta) * r_u;
									if (dtheta == 0) {
										x_vu = fabs(r_u - r_v);
									} else {
										x_vu =
												(1 / zeta)
														* acosh(
																(cosh(
																		zeta
																				* r_v)
																		* cosh(
																				zeta
																						* r_u))
																		- (sinh(
																				zeta
																						* r_v)
																				* sinh(
																						zeta
																								* r_u)
																				* cos(
																						dtheta)));
									}
									R_init = u->get_R();

								} else {

									double r_v = v->getInitRadius();
									double r_u = Beta * u->getInitRadius()
											+ (1 - Beta) * r_v;
									if (dtheta == 0) {
										x_vu = fabs(r_u - r_v);
									} else {
										x_vu =
												(1 / zeta)
														* acosh(
																(cosh(
																		zeta
																				* r_v)
																		* cosh(
																				zeta
																						* r_u))
																		- (sinh(
																				zeta
																						* r_v)
																				* sinh(
																						zeta
																								* r_u)
																				* cos(
																						dtheta)));
									}
									R_init = v->get_R();

								}

								double P_uv =
										1
												/ (1
														+ exp(
																(zeta / (2 * T))
																		* (x_vu
																				- R_init)));

								//If flags are set on both nodes we ignore the links between them (see line 180).
								if ((u->getFlag() == 1)
										&& (v->getFlag() == 1)) {
									log_L = log_L + log(1 - P_uv);
									continue;
								}

								if (u->checkLink(v)) {
									log_L = log_L + log(P_uv);
								} else {
									log_L = log_L + log(1 - P_uv);
								}
							}
						}

						if (log_L >= maxlog_L) {
							u->setAngle(theta_u);
							maxlog_L = log_L;
						}
						theta = theta + step;
					}
				}	    //Endfor (pick next node).
			}	    //Endfor round.
		}	    //Endif for correction step.
	}	    //Endfor (pick next node degree).

//Output the node coordinates in the form: node_id angular_coordinate radial_coordinate.
	auto nodes2 = model->getNodeList();
	for (std::size_t j = 0; j < nodes2.size(); j++) {
		u = nodes2[j];
		nodesCoord[u->getId()] = std::make_pair(u->getAngle(), u->getRadius());
		logger(logDebug,
				std::to_string(u->getId()) + " " + std::to_string(u->getAngle())
						+ " " + std::to_string(u->getRadius()))
	}
}

template<typename Network, typename EdgeRndIt,
		typename ScoreRndIt, typename EdgeRndOutIt> void UHYPPredictor<
		Network, EdgeRndIt, ScoreRndIt,
		EdgeRndOutIt>::predict(EdgeRndIt begin,
		EdgeRndIt end, ScoreRndIt oscores) {
	logger(logDebug, "Predicting links...")
#ifdef LINKPRED_WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		auto coordsi = nodesCoord.at(Network::start(*it));
		auto coordsj = nodesCoord.at(Network::end(*it));

		double dtheta = std::abs(coordsi.first - coordsj.first);
		if (dtheta > MathPI) {
			dtheta = 2 * MathPI - dtheta;
		}
		double ri = coordsi.second;
		double rj = coordsj.second;
		double dist = (1 / zeta)
				* std::acosh(
						(std::cosh(zeta * ri) * std::cosh(zeta * rj))
								- (std::sinh(zeta * ri) * std::sinh(zeta * rj)
										* std::cos(dtheta)));
		*(oscores + (it - begin)) = 1 / (1 + std::exp((zeta / (2 * T)) * dist));
	}
	logger(logDebug, "Done")
}

#define UHYPPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef UHYPPREDICTOR_CPP


} /* namespace LinkPred */
