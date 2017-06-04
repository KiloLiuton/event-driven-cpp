#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <string>
#include <sstream>

#include "pcg_random.hpp"
#include "topology.hpp"
#include "lattice.hpp"

#define NON_DETERMINISTIC_SEED true

// TODO:
// - monitor file changes with python script for real time plotting
// - write a better README.md using the markdown language

int main(int argc, char *argv[]) {

	// define lattice parameters:
	// any changes regarding topology should be done by creating a new lattice instance.
	const int SIZE = 201;
	const int K = 50;
	const double REWIRE_PROB = 0.0; // TODO : this is not implemented yet! only create regular rings
	double couplingStrength = 2.0;

	// set simulation parameters
	const size_t ITERS = 20000; // maximum amount of iterations in case system takes too long to relax
	const size_t TRIALS = 300; // number of independent runs for each 'couplingStrength' value

	// paths to data storage folders
	std::string rvsaData ("../rvsaData/");
	std::string relaxationData ("../relaxationData/");

	// create relaxation&rvsa filenames. If either exists, append a '*' to its name
	std::ostringstream oss;
	oss << "relaxation-" << "N=" << SIZE << "k=" << K << "p=" << REWIRE_PROB
		<< "a=" << couplingStrength << "TRIALS=" << TRIALS << ".txt";
	std::string relaxationFilename = oss.str();
	while(std::ifstream(relaxationData + relaxationFilename)) {
		relaxationFilename = relaxationFilename.substr(0, relaxationFilename.size()-4) + "+.txt";
	}
	oss.str("");
	oss << "rvsa-" << "N=" << SIZE << "k=" << K << "p=" << REWIRE_PROB
		<< "TRIALS=" << TRIALS << ".txt";
	std::string rvsaFilename = oss.str();
	while(std::ifstream(rvsaData + rvsaFilename)) {
		rvsaFilename = rvsaFilename.substr(0, rvsaFilename.size()-4) + "+.txt";
	}

	// open the new files for writing
	std::ofstream relaxationFile (relaxationData + relaxationFilename);
	std::ofstream rvsaFile (rvsaData + rvsaFilename);
	if(!relaxationFile.is_open())
		throw std::runtime_error("failed to open relaxation file. Make sure 'relaxationData' folder exists.");
	if(!rvsaFile.is_open())
		throw std::runtime_error("failed to open rvsa file. Make sure 'rvsaData' folder exists.");

	// seed rng
	pcg64 rng(42u, 54u);
	if(NON_DETERMINISTIC_SEED) rng.seed(pcg_extras::seed_seq_from<std::random_device>());

	// CREATE LATTICE INSTANCE
	Lattice lattice(SIZE, K, REWIRE_PROB, couplingStrength, rng);

	// write relaxatoin header
	relaxationFile << "# data used to determine relaxation period.\n"
		           << "# dt\tr\n";

	// relaxation run
	for(size_t i = 0; i < ITERS; ++i) {
		double dt = lattice.step();
		double r = lattice.getOrderParameter();
		relaxationFile << dt << "\t" << r << std::endl;
	}

	// FIXME: detect relaxation and set relaxationPeriod and pointsAfterRelaxation to sensible values
	size_t relaxationPeriod = 0;
	size_t pointsAfterRelaxation = 666666666666;

	// r vs a run
	int numPoints = 20;
	double initialCoupling = 1.3, finalCoupling = 3.6;
	std::vector<double> aRange (numPoints);
	for(int i = 0; i < numPoints; ++i) {
		aRange[i] = initialCoupling + i * (finalCoupling - initialCoupling) / numPoints;
	}

	// write rvsa header
	rvsaFile << "# TRIALS=" << TRIALS << "\trelaxationPeriod=" << relaxationPeriod
	         << "\tpointsAfterRelaxation=" << pointsAfterRelaxation << std::endl
			 << "# a" << "\t<<r>>" << "\tX=<<r2>>-<<r>>2\tX'=<<r>2>-<<r>>2\n";

	// begin coupling strength set
	for(size_t a = 0; a < aRange.size(); ++a) {
		lattice.setCouplingStrength(aRange[a]);
		double rAvgSum = 0;
		double r2AvgSum = 0;
		double rAvg2Sum = 0;
		for(size_t j = 0; j < TRIALS; ++j) { // run TRIALS trials for each coupling strength
			double rSum = 0;
			double r2Sum = 0;
			double dtSum = 0;
			// each trial consists of 'ITERS' steps at most, and the first 'relaxationPeriod' steps are discarded
			for(size_t i = 0; i < ITERS && i < relaxationPeriod + pointsAfterRelaxation; ++i) {
				double dt = lattice.step();
				if(i >= relaxationPeriod) {
					double r = lattice.getOrderParameter();
					rSum += r*dt;
					r2Sum += r*r*dt;
					dtSum += dt;
				}
			}
			double rTrialAvg = rSum / dtSum;
			rAvgSum += rTrialAvg;
			r2AvgSum += r2Sum / dtSum;
			rAvg2Sum += rTrialAvg*rTrialAvg;
		}
		double rAvgAvg = rAvgSum / TRIALS;
		double r2AvgAvg = r2AvgSum / TRIALS;
		double rAvg2Avg = rAvg2Sum / TRIALS;
		double X = r2AvgAvg - rAvgAvg*rAvgAvg;
		double Xnew = rAvg2Avg - rAvgAvg*rAvgAvg;
		std::cout << aRange[a] << " finished\t" << "[" << a+1 << "/" << numPoints << "]\n";
		rvsaFile << std::fixed << std::setprecision(10)
		         << aRange[a] << "\t" << rAvgAvg << "\t" << X << "\t" << Xnew << std::endl;
	}

	return 0;
}
