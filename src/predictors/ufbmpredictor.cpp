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

#include <linkpred/predictors/ufbmpredictor.hpp>
#include <cmath>
#include <omp.h>
#include "linkpred/utils/utilities.hpp"

namespace LinkPred {

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%% train: adjacency matrix of the observed network.
//%%% n: number of nodes of the observed network.
//%%% i: the ith node in the network.
//%%% cl(output): a clique conludes node i.
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//function  cl  = clique_find(train,n,i)
//    a=train(i,:);
//    b=find(a);
//    c=train(b',:);
//    ct=1;
//    cl(ct)=i;
//    while (sum(a)>0)
//        d=c*a';
//        md=max(d);
//        p=find(d==md);
//        i=b(p(1));
//        a=c(p(1),:).*a;
//        b=find(a);
//        c=train(b',:);
//        ct=ct+1;
//        cl(ct)=i;
//    end
//end

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> std::vector<
		std::size_t> UFBMPredictor<NetworkT, EdgesRandomIteratorT,
		ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::cliqueFind(
		GFMatrix const & train, std::size_t i) {
	Vec a = train.getRow(i);
//	std::cout << "a: " << std::endl;
//	a.print();
	auto b = a.findNz();
//	std::cout << "b: " << std::endl;
//	for (auto it = b.begin(); it != b.end(); ++it) {
//		std::cout << *it << std::endl;
//	}
	GFMatrix c = train.getRows(b);
//	std::cout << "c: " << std::endl;
//	c.print();
	std::vector < std::size_t > cl;
	cl.push_back(i);
	while (a.sum() > 0) {
		Vec d = c * a;
//		std::cout << "d: " << std::endl;
//		d.print();
		double md = d.max();
//		std::cout << "md: " << md << std::endl;
		auto p = d.find(md);
//		std::cout << "p: " << std::endl;
//		for (auto it = p.begin(); it != p.end(); ++it) {
//			std::cout << *it << std::endl;
//		}
		i = b[p[0]];
//		std::cout << "i: " << i << std::endl;
//		c.getRow(p[0]).print();
		a = c.getRow(p[0]) * a; // * stands here for element-wise product
//		std::cout << "a: " << std::endl;
//		a.print();
		b = a.findNz();
		c = train.getRows(b);
		cl.push_back(i);
	}
	return cl;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%% s1: the node vector of the current netowork.
//%%% relat: adjacency matrix of the current network.
//%%% n: nubmer of nodes.
//%%% s1(output): an obtained node vector which denots a clique.
//%%% s2(output): a node vertor which denotes the remainder of the network.
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//function [s1,s2]=block_g(s1,relat,n)
//    dg=sum(relat);
//    [sdg lx]=sort(dg);
//    k=lx(n);
//    cl = clique_find(relat,n,k);
//    s2=s1;
//    s2(cl)=0;
//    s1=s1-s2;
//end

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> std::pair<
		std::vector<std::size_t>, std::vector<std::size_t>> UFBMPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::blockG(std::vector<std::size_t> const & s1,
		GFMatrix const & relat) {
	Vec dg = relat.sumCols();
	std::size_t k = dg.ilmax();
	auto cl = cliqueFind(relat, k);
	auto s2 = s1;
	for (auto it = cl.begin(); it != cl.end(); ++it) {
		s2[*it] = 0;
	}
	std::vector < std::size_t > ns1;
	ns1.resize(s2.size());
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (std::size_t i = 0; i < s2.size(); i++) {
		ns1[i] = s1[i] - s2[i];
	}
	return std::make_pair(ns1, s2);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%% train: adjacency matrix of the observed network.
//%%% n: number of nodes of the observed network.
//%%% Clus: obtained communities.
//%%% k: number of communities.
//%%% ts: the index of the first special community in Clus
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//function [clus k ts]=division_g(relat,n)
//    %%%%%%%%%%%%%%%%%randomly divide the network into two parts%%%%%%%%%%%%%%%%%%
//    k=2;
//    x=randperm(n);
//    t=floor(n/k);
//    for (i=1:k)
//        gp(i,:)=zeros(1,n);
//        if (i<k)
//            s=x(t*(i-1)+1:t*i);
//        else
//            s=x(t*(i-1)+1:n);
//        end
//        gp(i,s)=1;
//    end
//    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//    clus=zeros(1,n);
//    j=0;
//    for (i=1:k)
//        relat1=(gp(i,:)'*gp(i,:)).*relat;
//        l=sum(sum(relat1))/2;
//        relat2=relat1;
//        while (l>0)
//             [s1,gp(i,:)]=block_g(gp(i,:),relat2,n);
//             j=j+1;
//             clus(j,:)=s1;
//             relat2=(gp(i,:)'*gp(i,:)).*relat1;
//             n1=sum(gp(i,:));
//             l=sum(sum(relat2))/2;
//        end
//        j=j+1;
//        clus(j,:)=gp(i,:);
//        if (i==1)
//            ts=j;
//        end
//    end
//    k=j;
//end

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> std::pair<
		GFMatrix, std::size_t> UFBMPredictor<NetworkT, EdgesRandomIteratorT,
		ScoresRandomIteratorT, EdgesRandomOutputIteratorT>::divisionG(
		GFMatrix const & train) {
	// randomly divide the network into two part
	std::size_t ts;
	std::size_t n = train.getN();
	std::size_t t = std::floor((double) n / 2);
	auto x = Utilities::getRndPerm(n, rng.getInt());
//	Utils::print(x, "x");
//	std::vector<std::size_t> x = { 4, 33, 15, 21, 1, 26, 29, 32, 30, 7, 18, 25, 0, 17,
//			24, 10, 2, 19, 5, 13, 22, 23, 6, 31, 3, 8, 16, 20, 11, 12, 27, 9,
//			28, 14 };
//	std::vector<std::size_t> x = { 38, 10, 47, 41, 33, 2, 63, 26, 43, 14, 39, 29, 8, 23,
//			27, 7, 55, 42, 60, 4, 46, 53, 17, 28, 11, 18, 48, 34, 54, 61, 0, 58,
//			57, 12, 22, 59, 25, 1, 44, 45, 15, 36, 6, 56, 9, 19, 21, 13, 50, 35,
//			40, 37, 31, 30, 20, 51, 49, 62, 52, 3, 16, 24, 5, 32 };
//	Utils::print(x, "x");

	GFMatrix gp(2, Utilities::int_cast(n), true);
//	gp.setName("gp");
	for (std::size_t s = 0; s < t; s++) {
		gp.set(0, x[s], 1);
	}
	for (std::size_t s = t; s < n; s++) {
		gp.set(1, x[s], 1);
	}

//	gp.print("gp");
	std::vector < std::vector < std::size_t >> clusVec;
	// ********* i = 0 *************
	{

		auto relat1 = GFMatrix::elemMult(
				GFMatrix::mult(gp.getRow(0), gp.getRow(0)), train);
//		GFMatrix relat1(n, n);
//		for (std::size_t i = 0; i < n; i++) {
//			for (std::size_t j = 0; j < n; j++) {
//				relat1.set(i, j, gp.get(0, i) * gp.get(0, j) * train.get(i, j));
//			}
//		}
//		relat1.print("relat1");

		double l = relat1.sum() / 2;

//		std::cout << "l: " << l << std::endl;
		GFMatrix relat2(relat1);
//		gp.setName("relat2");

		while (l > 0) {

			auto row = gp.getRow(0);

			std::vector < std::size_t > vec;
			for (std::size_t i = 0; i < n; i++) {
				vec.push_back(row[i]);
			}

			auto res = blockG(vec, relat2);

//			Utils::print(res.first, "s1");
//			Utils::print(res.second, "s2");
			for (std::size_t j = 0; j < n; j++) {
				gp.set(0, j, res.second[j]);
			}

			clusVec.push_back(res.first);

//			std::cout << "clus:" << std::endl;
//			for (auto it1 = clus.begin(); it1 != clus.end(); ++it1) {
//				for (auto it2 = it1->begin(); it2 != it1->end(); ++it2) {
//					std::cout << *it2 << "\t";
//				}
//				std::cout << std::endl;
//			}

//             relat2=(gp(i,:)'*gp(i,:)).*relat1;
			relat2 = GFMatrix::elemMult(
					GFMatrix::mult(gp.getRow(0), gp.getRow(0)), relat1);
//			for (std::size_t i = 0; i < n; i++) {
//				for (std::size_t j = 0; j < n; j++) {
//					relat2.set(i, j,
//							gp.get(0, i) * gp.get(0, j) * relat1.get(i, j));
//				}
//			}
//			gp.print("gp");
			//             l=sum(sum(relat2))/2;

			l = relat2.sum() / 2;
//			std::cout << "i = 0, l: " << l << std::endl;
		}

		auto row = gp.getRow(0);
		std::vector < std::size_t > vec;
		for (std::size_t i = 0; i < n; i++) {
			vec.push_back(row[i]);
		}

		clusVec.push_back(vec);
		ts = clusVec.size() - 1;
	}

	// ********* i = 1 *************
	{
		auto relat1 = GFMatrix::elemMult(
				GFMatrix::mult(gp.getRow(1), gp.getRow(1)), train);
//		GFMatrix relat1(n, n);
//		for (std::size_t i = 0; i < n; i++) {
//			for (std::size_t j = 0; j < n; j++) {
//				relat1.set(i, j, gp.get(1, i) * gp.get(1, j) * train.get(i, j));
//			}
//		}
//		relat1.print("relat1");

		double l = relat1.sum() / 2;
//		std::cout << "l: " << l << std::endl;

		GFMatrix relat2(relat1);
//		gp.setName("relat2");

		while (l > 0) {

			auto row = gp.getRow(1);
			std::vector < std::size_t > vec;
			for (std::size_t i = 0; i < n; i++) {
				vec.push_back(row[i]);
			}

			auto res = blockG(vec, relat2);

//			Utils::print(res.first, "s1");
//			Utils::print(res.second, "s2");
			for (std::size_t j = 0; j < n; j++) {
				gp.set(1, j, res.second[j]);
			}

			clusVec.push_back(res.first);

//			std::cout << "clus:" << std::endl;
//			for (auto it1 = clus.begin(); it1 != clus.end(); ++it1) {
//				for (auto it2 = it1->begin(); it2 != it1->end(); ++it2) {
//					std::cout << *it2 << "\t";
//				}
//				std::cout << std::endl;
//			}

			//             relat2=(gp(i,:)'*gp(i,:)).*relat1;
			relat2 = GFMatrix::elemMult(
					GFMatrix::mult(gp.getRow(1), gp.getRow(1)), relat1);
//			for (std::size_t i = 0; i < n; i++) {
//				for (std::size_t j = 0; j < n; j++) {
//					relat2.set(i, j,
//							gp.get(1, i) * gp.get(1, j) * relat1.get(i, j));
//				}
//			}
			//             l=sum(sum(relat2))/2;
//			gp.print("gp");

			l = relat2.sum() / 2;

//			std::cout << "i = 1, l: " << l << std::endl;
		}

		auto row = gp.getRow(1);

		std::vector < std::size_t > vec;
		for (std::size_t i = 0; i < n; i++) {
			vec.push_back(row[i]);
		}

		clusVec.push_back(vec);
	}

	GFMatrix clus(Utilities::int_cast(clusVec.size()), Utilities::int_cast(n));

	for (std::size_t i = 0; i < clusVec.size(); i++) {
		for (std::size_t j = 0; j < n; j++) {
			clus.set(i, j, clusVec[i][j]);
		}
	}
	return std::make_pair(clus, ts);
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UFBMPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::init() {
//create adj matrix
	logger(logDebug, "Initializing solver...")
#ifdef WITH_OPENMP
	if (parallel) {
//		throw std::runtime_error(
//				"Parallelism is not yet properly supported in FBM");
		std::cerr << "# Warning: Parallelism is not yet properly supported in FBM. It will be disabled\n";
		parallel = false;
	}
#endif
	auto n = net->getNbNodes();
	adj = std::make_shared < GFMatrix
			> (Utilities::int_cast(n), Utilities::int_cast(n), true); // Create and init to 0

//#pragma omp parallel for
	for (auto it = net->edgesBegin(); it < net->edgesEnd(); ++it) {
		adj->set(NetworkT::start(*it), NetworkT::end(*it), 1);
		adj->set(NetworkT::end(*it), NetworkT::start(*it), 1);
	}
	logger(logDebug, "Done")
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%% train: adjacency matrix of the observed network.
//%%% n: number of nodes of the observed network.
//%%% W: a symmetrical matrix containing the link probabilities for all non-adjacent
//%%% links in the observed network.
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//function  W  = FBM( train,n)
//    rs=zeros(n,n);
//    tt=50; % sampling parameter
//    for (ii=1:tt)
//         [clus, k, ts] = division_g(train,n);
//         dg=sum(clus');
//         c2=clus*train*clus';
//         diag1=c2.*eye(k)/2;
//         c2=c2-c2.*eye(k);
//         c3=dg'*dg;
//         c3=c3-c3 .* eye(k);
//         diag2=dg'*(dg-1)./2.*eye(k);
//         diag=diag2.*2-diag1;
//         c3=c3+c2;
//         c3=c3+diag;
//         diag=diag2;
//         diag(ts,ts)=0;
//         diag(k,k)=0;
//         c2=c2+diag;
//         r_block=c2./c3;
//         rs=rs+clus'*r_block*clus*k;
//    end
//    dig = rs .* eye(size(rs));
//    rs = (rs - dig)./tt;
//    W = rs-rs.*train;
//end
template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UFBMPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::learn() {
	logger(logDebug, "Learning...")
#ifdef WITH_OPENMP
	if (parallel) {
		throw std::runtime_error(
				"Parallelism is not yet properly supported in FBM");
	}
#endif
	fbm (adj);
	logger(logDebug, "Done")
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UFBMPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::fbm(std::shared_ptr<GFMatrix> adj) {

//	adj->print("Adj");
	std::size_t n = adj->getN();
	GFMatrix rs(Utilities::int_cast(n), Utilities::int_cast(n), true);

//	maxIter = 1;
	for (std::size_t iter = 0; iter < maxIter; iter++) {
//		std::cout << "iter: " << iter << std::endl;
		auto res = divisionG(*adj);

//		res.first.print("clus");
		auto dg = res.first.sumRows();
//		dg.print("dg");

		auto c2 = GFMatrix::mult(res.first * (*adj), res.first, false, true);
//		c2.print("c2");
//		exit(1);

		Vec diag1 = c2.diag() / 2;
//		diag1.print("diag1");

		c2.setDiag(0);
//		c2.print("c2");

		auto c3 = GFMatrix::mult(dg, dg);
//		c3.print("c3");

		c3.setDiag(0);
//		c3.print("c3");

		Vec diag2(dg.size());

		for (int i = 0; i < dg.size(); i++) {
			diag2[i] = dg[i] * (dg[i] - 1) / 2; // >=0
		}
//		diag2.print("diag2");

		Vec diag = 2 * diag2 - diag1;
//		diag.print("diag");

		c3 = c3 + c2;
//		c3.print("c3");
		c3 = c3 + diag;
//		c3.print("c3");

		diag = diag2;
//		diag.print("diag");

		diag[res.second] = 0;
//		diag.print("diag");

		diag[res.first.getM() - 1] = 0;
//		diag.print("diag");

		c2 = c2 + diag;
//		c2.print("c2");

//		c3.print("c3");
		auto rBlock = c2 / c3; // element-wise division
		rBlock.removeNaN();
//		rBlock.print("rBlock");

		rs = rs
				+ res.first.getM()
						* GFMatrix::mult(res.first, rBlock, true, false)
						* res.first;
	}

	rs.setDiag(0);
	rs = (1.0 / maxIter) * rs;
	rs = rs - GFMatrix::elemMult(rs, *adj);

//	rs.print("rs");

	for (auto it = net->nonEdgesBegin(); it != net->nonEdgesEnd(); ++it) {
		auto i = NetworkT::start(*it);
		auto j = NetworkT::end(*it);
		double w = adj->get(i, j);
		if ((w != 0) && (w != 1)) {
			std::cout << i << "\t" << j << "\t" << adj->get(i, j) << "\t"
					<< net->getNbNodes() << std::endl;
			throw std::runtime_error("Invalid value in the adjacency matrix");
		}
		if (w == 1) {
			throw std::runtime_error(
					"Corrupted adjacency matrix. This entry should be 0");
		}
		scores[*it] = rs.get(i, j);
	}
}

template<typename NetworkT, typename EdgesRandomIteratorT,
		typename ScoresRandomIteratorT, typename EdgesRandomOutputIteratorT> void UFBMPredictor<
		NetworkT, EdgesRandomIteratorT, ScoresRandomIteratorT,
		EdgesRandomOutputIteratorT>::predict(EdgesRandomIteratorT begin,
		EdgesRandomIteratorT end, ScoresRandomIteratorT oscores) {
	logger(logDebug, "Predicting links...")
#ifdef WITH_OPENMP
#pragma omp parallel for if (parallel)
#endif
	for (auto it = begin; it < end; ++it) {
		*(oscores + (it - begin)) = scores.at(*it);
	}
	logger(logDebug, "Done")
}

#define UFBMPREDICTOR_CPP
#include "linkpred/instantiations.hpp"
#undef UFBMPREDICTOR_CPP

}
/* namespace LinkPred */
