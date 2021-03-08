#include <linkpred.hpp>
#include <iostream>

using namespace LinkPred::Simp;

int main() {
#ifdef LINKPRED_WITH_MLPACK
	int k = 10;
	Predictor p;
	p.loadnet("Zakarays_Karate_Club.edges");
	{
		auto esv = p.predTopECL(k, "DPW", "FFN");
		std::cout << "Top " << k << " using DPW-FFN\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "DPW", "LSVM");
		std::cout << "Top " << k << " using DPW-LSVM\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "DPW", "LGR");
		std::cout << "Top " << k << " using DPW-LGR\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "DPW", "NVB");
		std::cout << "Top " << k << " using DPW-NVB\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}

	{
		auto esv = p.predTopECL(k, "HMSM", "FFN");
		std::cout << "Top " << k << " using HMSM-FFN\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "HMSM", "LSVM");
		std::cout << "Top " << k << " using HMSM-LSVM\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "HMSM", "LGR");
		std::cout << "Top " << k << " using HMSM-LGR\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "HMSM", "NVB");
		std::cout << "Top " << k << " using HMSM-NVB\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}

	{
		auto esv = p.predTopECL(k, "LVS", "FFN");
		std::cout << "Top " << k << " using LVS-FFN\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "LVS", "LSVM");
		std::cout << "Top " << k << " using LVS-LSVM\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "LVS", "LGR");
		std::cout << "Top " << k << " using LVS-LGR\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "LVS", "NVB");
		std::cout << "Top " << k << " using LVS-NVB\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}

	{
		auto esv = p.predTopECL(k, "LIN", "FFN");
		std::cout << "Top " << k << " using LIN-FFN\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "LIN", "LSVM");
		std::cout << "Top " << k << " using LIN-LSVM\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "LIN", "LGR");
		std::cout << "Top " << k << " using LIN-LGR\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "LIN", "NVB");
		std::cout << "Top " << k << " using LIN-NVB\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}

	{
		auto esv = p.predTopECL(k, "N2V", "FFN");
		std::cout << "Top " << k << " using N2V-FFN\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "N2V", "LSVM");
		std::cout << "Top " << k << " using N2V-LSVM\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "N2V", "LGR");
		std::cout << "Top " << k << " using N2V-LGR\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
	{
		auto esv = p.predTopECL(k, "N2V", "NVB");
		std::cout << "Top " << k << " using N2V-NVB\n";
		for (auto it = esv.begin(); it != esv.end(); ++it) {
			std::cout << it->i << "\t" << it->j << "\t" << it->score
					<< std::endl;
		}
	}
#endif

	return 0;
}
