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
	const int SIZE = 400;
	const int K = 60;
	const double REWIRE_PROB = 0.0; // TODO : this is not implemented yet! only create regular rings
	double couplingStrength = 1.3;

	// set simulation parameters
	const size_t MAX_ITERS = 1000; // maximum amount of iterations in case system takes too long to relax
	const size_t TRIALS = 200; // number of independent runs for each 'couplingStrength' value

	// paths to data storage folders
	std::string rvsaData ("rvsaData/");
	std::string relaxationData ("relaxationData/");

	// create relaxation&rvsa filenames. If either exists, append a '+' to its name
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
	Lattice simulation(SIZE, K, REWIRE_PROB, couplingStrength, rng);
	simulation.printTopology();

	// relaxation run
	// write relaxation header
	relaxationFile << "# data used to determine relaxation period.\n"
		           << "# dt\tr\tN0\tN1\n";

	// FIXME: This function might not be the best solution to detecting relaxation.
	//	   This function gets the average of the order parameter for a block of 'trail' events. The next block
	//		   is obtained by shifting the block one event (discard the first event and include the next).
	//		   If the average of the new block is less than 'threshold' of the average of the first block,
	//		   increment a count value. If this count value reaches a size equal to trail, break the loop and
	//		   return the current value of steps.
	// TODO: simulation run currently ignores any values found by the relaxation period step. Implement
	//       this after the relaxation period has been fixed.
	//
	// size_t relaxationRun(trail, threshold, MAX_ITERS, outputfile) threshold>1.0 ensures it runs for MAX_ITERS
	//     the order parameter, step duration and populations are recorded in outputfile for every step.
	size_t relaxationPeriod = simulation.relaxationRun(100, 2.0, MAX_ITERS, relaxationFile);
	size_t pointsAfterRelaxation = 666666666666;

	// r vs a run
	// start by creating a vector with the coupling strength values
	int numPoints = 20;
	double initialCoupling = 1.3, finalCoupling = 3.6;
	std::vector<double> aRange (numPoints);
	for(int i=0; i<numPoints; ++i) aRange[i] = initialCoupling + i*(finalCoupling - initialCoupling)/numPoints;

	// write rvsa header (order parameter r vs coupling strength a)
	rvsaFile << "# TRIALS=" << TRIALS << "\trelaxationPeriod=" << relaxationPeriod
	         << "\tpointsAfterRelaxation=" << pointsAfterRelaxation << std::endl
			 << "# a" << "\t<<r>>" << "\tX=<<r2>>-<<r>>2\tX'=<<r>2>-<<r>>2\n";

	// TODO: isolate trial-run into it's own function in lattice.cpp (to simplify the nested loops)
	//
	// outer loop: set coupling strength
	//     middle loop: start a trial
	//         inner loop: run a single trial for MAX_ITERS steps at most
	for(size_t a = 0; a < aRange.size(); ++a) {
		simulation.setCouplingStrength(aRange[a]);
		double rAvgSum = 0;
		double r2AvgSum = 0;
		double rAvg2Sum = 0;
		for(size_t j = 0; j < TRIALS; ++j) { // run TRIALS trials for each coupling strength
			double rSum = 0;
			double r2Sum = 0;
			double dtSum = 0;

			// each trial consist of 'MAX_ITERS' steps at most, and the first 'relaxationPeriod' steps are discarded
			for(size_t i = 0; i < MAX_ITERS; ++i) {
				double dt = simulation.step();
				double r = simulation.getOrderParameter();
				rSum += r*dt;
				r2Sum += r*r*dt;
				dtSum += dt;
			}
			double rAvg = rSum / dtSum;
			double r2Avg = r2Sum / dtSum;

			rAvgSum += rAvg;
			r2AvgSum += r2Avg;
			rAvg2Sum += rAvg*rAvg;
		}
		double rAvgAvg = rAvgSum / TRIALS;        // <<r>>
		double r2AvgAvg = r2AvgSum / TRIALS;      // <<r^2>>
		double rAvg2Avg = rAvg2Sum / TRIALS;      // <<r>^2>
		double X = r2AvgAvg - rAvgAvg*rAvgAvg;    // <<r^2>> - <<r>>^2
		double Xnew = rAvg2Avg - rAvgAvg*rAvgAvg; // <<r>^2> - <<r>>^2

		std::cout << aRange[a] << " finished\t" << "[" << a+1 << "/" << numPoints << "]\n"; // track progress
		rvsaFile << std::fixed << std::setprecision(12)
		         << aRange[a] << "\t" << rAvgAvg << "\t" << X << "\t" << Xnew << std::endl;
	}

	return 0;
}
