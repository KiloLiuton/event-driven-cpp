#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <sstream>

#include "pcg_random.hpp"
#include "topology.h"
#include "lattice.h"

#define non_deterministic_seed false

int main(int argc, char *argv[]) {

	// define lattice parameters:
	// any changes regarding topology should be done by creating a new lattice instance.
	const int SIZE = 200;
	const int K = 20;
	const double REWIRE_PROBABILITY = 0.0; // TODO : this is not implemented yet! only regular rings are created
	double couplingStrength = 1.6;

	// set simulation parameters
	const size_t ITERS = 200000; // maximum amount of iterations in case system takes too long to relax
	const size_t TRIALS = 1; // number of independent runs for each 'couplingStrength' value

	// prepare output filenames
	std::ostringstream oss;
	oss << "N=" << SIZE << "k=" << K << "p=" << REWIRE_PROBABILITY << "a=" << couplingStrength << "ITERS=" << ITERS << "TRIALS=" << TRIALS << ".txt";
	std::string relaxationFilename = oss.str();

	oss << "N=" << SIZE << "k=" << K << "p=" << REWIRE_PROBABILITY << "ITERS=" << ITERS << "TRIALS=" << TRIALS << ".txt";
	std::string rvsaFilename = oss.str();

	// open files and write headers
	std::string relaxationDataPath ("relaxationData/"); // folder name for relaxation data
	std::string rvsaDataPath ("rvsaData/"); // folder name for rvsa data
	std::ofstream relaxationFile (relaxationDataPath + relaxationFilename);
	if(!relaxationFile.is_open()) throw std::runtime_error("failed to open relaxation file. Make sure 'relaxationData' folder exists.");
	std::ofstream rvsaFile (rvsaDataPath + rvsaFilename);
	if(!rvsaFile.is_open()) throw std::runtime_error("failed to open rvsa file. Make sure 'rvsaData' folder exists.");

	relaxationFile << "# recording of 'r' for relaxation checking. Using 'a=" << couplingStrength << "'\n";
	relaxationFile << "# dt\tr\n";
	rvsaFile << "# average 'r' and 'X' values after relaxation period for each value of 'a'\n";
	rvsaFile << "# a\t<<r>>\tX=<<r2>>-<<r>>2\n";

	// seed rng
	pcg64 rng(42u, 54u);
	if(non_deterministic_seed) rng.seed(pcg_extras::seed_seq_from<std::random_device>());

	// create a lattice instance
	Lattice lattice(SIZE, K, REWIRE_PROBABILITY, couplingStrength, rng);

	// print topology and initial condition FIXME: (may print incorrectly if SIZE is too big)
	/*
	lattice.printTopology();
	std::cout << "\nINITIAL CONDITION\n";
	lattice.print();
	*/

	// relaxation run
	for(size_t i = 0; i < ITERS; ++i) {
		double dt = lattice.step();
		double r = lattice.getOrderParameter();
		relaxationFile << dt << "\t" << r << std::endl;
	}

	// FIXME: detect relaxation and set relaxationPeriod and relaxedData to sensible values
	size_t relaxationPeriod = 0;
	size_t relaxedData = 666666666666;

	// r vs a run
	double rAvgSum = 0;
	double r2AvgSum = 0;
	for(size_t j = 0; j < TRIALS; ++j) {
		double rSum = 0;
		double r2Sum = 0;
		double dtSum = 0;
		for(size_t i = 0; i < ITERS && i < relaxationPeriod + relaxedData; ++i) {
			double dt = lattice.step();
			if(i >= relaxationPeriod) {
				double r = lattice.getOrderParameter();
				double r2 = r*r;
				rSum += r*dt;
				r2Sum += r2*dt;
				dtSum += dt;
			}
		}
		rAvgSum += rSum/dtSum;
		r2AvgSum += r2Sum/dtSum;
	}
	double rAvgAvg = rAvgSum/TRIALS;
	double r2AvgAvg = r2AvgSum/TRIALS;
	double X = r2AvgAvg - rAvgAvg*rAvgAvg;
	std::cout << "<r>2 = " << rAvgAvg*rAvgAvg << "\t<r2> = " << r2AvgAvg << std::endl;
	rvsaFile << couplingStrength << "\t" << rAvgAvg << "\t" << X << std::endl;


	return 0;
}
