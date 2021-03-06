#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <string>
#include <sstream>

#include "pcg_random.hpp"
#include "topology.hpp"
#include "lattice.hpp"

static int LATTICE_SIZE = 801;
static int NUMBER_OF_FORWARD_NEIGHBORS = 50;
static float REWIRE_PROBABILITY = 0.075;
static int MAXIMUM_ITERATIONS = 25000;
static int NUMBER_OF_TRIALS = 600;
static int NON_DETERMINISTIC_SEED = 1;
static float RELAXATION_COUPLING = 2.59;
static int RELAXATION_BLOCK_SIZE = 100;
static float RELAXATION_THRESHOLD = 0.005;
static int TIMES_TO_RESET = 5;


// TODO:
// - monitor file changes with python script for real time plotting
// - write a better README.md using the markdown language

int main(int argc, char *argv[]) {
	if(auto tmp = getenv("LATTICE_SIZE")) { LATTICE_SIZE = atoi(tmp); }
	if(auto tmp = getenv("NUMBER_OF_FORWARD_NEIGHBORS")) { NUMBER_OF_FORWARD_NEIGHBORS = atoi(tmp); }
	if(auto tmp = getenv("REWIRE_PROBABILITY")) { REWIRE_PROBABILITY = atof(tmp); }
	if(auto tmp = getenv("MAXIMUM_ITERATIONS")) { MAXIMUM_ITERATIONS = atoi(tmp); }
	if(auto tmp = getenv("NUMBER_OF_TRIALS")) { NUMBER_OF_TRIALS = atoi(tmp); }
	if(auto tmp = getenv("NON_DETERMINISTIC_SEED")) { NON_DETERMINISTIC_SEED = atoi(tmp); }
	if(auto tmp = getenv("RELAXATION_COUPLING")) { RELAXATION_COUPLING = atof(tmp); }
	if(auto tmp = getenv("RELAXATION_BLOCK_SIZE")) { RELAXATION_BLOCK_SIZE = atoi(tmp); }
	if(auto tmp = getenv("RELAXATION_THRESHOLD")) { RELAXATION_THRESHOLD = atof(tmp); }
	if(auto tmp = getenv("TIMES_TO_RESET")) { TIMES_TO_RESET = atoi(tmp); }

	// define lattice parameters:
	// any changes regarding topology should be done by creating a new lattice instance.
	const int SIZE = LATTICE_SIZE;
	const int K = NUMBER_OF_FORWARD_NEIGHBORS;
	const double REWIRE_PROB = REWIRE_PROBABILITY;
	double couplingStrength = RELAXATION_COUPLING;

	// set simulation parameters
	// maximum amount of iterations in case system takes too long to relax
	const size_t MAX_ITERS = MAXIMUM_ITERATIONS;
	// number of independent runs for each 'couplingStrength' value

	// paths to data storage folders
	std::string rvsaData ("rvsaData/");
	std::string relaxationData ("relaxationData/");

	// create relaxation&rvsa filenames. If either exists, append a '+' to its name
	std::ostringstream oss;
	oss << "relaxation-" << "N=" << SIZE << "k=" << K << "p=" << REWIRE_PROB
		<< "a=" << couplingStrength << "TRIALS=" << NUMBER_OF_TRIALS << "ITER=" << MAX_ITERS << ".txt";
	std::string relaxationFilename = oss.str();
	while(std::ifstream(relaxationData + relaxationFilename)) {
		relaxationFilename = relaxationFilename.substr(0, relaxationFilename.size()-4) + "+.txt";
	}
	oss.str("");
	oss << "rvsa-" << "N=" << SIZE << "k=" << K << "p=" << REWIRE_PROB
		<< "TRIALS=" << NUMBER_OF_TRIALS << "ITER=" << MAX_ITERS << ".txt";
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
	Lattice simulation(SIZE, K, REWIRE_PROB, false, couplingStrength, rng);
	//simulation.printTopology();

	// relaxation run
	// write relaxation header
	relaxationFile << "# data used to determine relaxation period.\n"
		           << "# dt\tr\tN0\tN1\n";

	// FIXME: This function might not be the best solution to detecting relaxation.
	//  This function gets the average of the order parameter for a block of 'trail' events. The next block
	// is obtained by shifting the block one event (discard the first event and include the next).
	//  If the average of the new block is less than 'threshold' of the average of the first block,
	// increment a count value. If this count value reaches a size equal to trail, break the loop and
	// return the current value of steps.

	// relaxationRun parameters:
	// blockSize - track the average OP in blockSize steps. Smalls values -> higher fluctuations
	size_t relaxationPeriod = simulation.relaxationRun(
			RELAXATION_BLOCK_SIZE,
			RELAXATION_THRESHOLD,
			2*MAX_ITERS,
			relaxationFile
			);
	size_t pointsAfterRelaxation = MAX_ITERS;
	std::cout << "Relaxation returned " << relaxationPeriod << " iterations for relaxation period.\n";
	std::cout << "Proceeding to burn " << relaxationPeriod << " steps and record " << pointsAfterRelaxation;
	std::cout << " for " << NUMBER_OF_TRIALS << " trials.\n";

	// r vs a run
	// start by creating a vector with the coupling strength values
	int numPoints = 20;
	double initialCoupling = 1.3, finalCoupling = 3.6;
	std::vector<double> aRange (numPoints);
	for(int i=0; i<numPoints; ++i) aRange[i] = initialCoupling + i*(finalCoupling - initialCoupling)/numPoints;

	// write rvsa header (order parameter r vs coupling strength a)
	rvsaFile << "# TRIALS=" << NUMBER_OF_TRIALS << "\trelaxationPeriod=" << relaxationPeriod
	         << "\tpointsAfterRelaxation=" << pointsAfterRelaxation << std::endl
			 << "# a" << "\t<<r>>" << "\tX=<<r2>>-<<r>>2\tX'=<<r>2>-<<r>>2\n";

	// TODO: isolate trial-run into it's own function in lattice.cpp (to simplify the nested loops)

	// outer loop: set coupling strength for trials
	//     middle loop: perform all trials
	//         inner loop: perform a single trial
	for(size_t a = 0; a < aRange.size(); ++a) {
		simulation.setCouplingStrength(aRange[a]);
		double rAvgSum = 0;
		double r2AvgSum = 0;
		double rAvg2Sum = 0;
		for(size_t j = 0; j < NUMBER_OF_TRIALS; ++j) { // run trials for each coupling strength
			double rSum = 0;
			double r2Sum = 0;
			double dtSum = 0;

			// discard the first 'relaxationPeriod' steps
			for (size_t i = 0; i < relaxationPeriod; ++i) simulation.step();

			// record data after relaxation period and for pointsAfterRelaxation
			// here we break up the loop in smaller chunks in order to refresh
			// the total rate and avoid numerical errors
			size_t chunkSize = pointsAfterRelaxation/TIMES_TO_RESET;
			for(size_t i = 0; i < TIMES_TO_RESET; ++i) {
				for (size_t j = 0; j < chunkSize; ++j) {
					double dt = simulation.step();
					double r = simulation.getOrderParameter();
					rSum += r*dt;
					r2Sum += r*r*dt;
					dtSum += dt;
				}
				simulation.resetTotalRate();
			}
			double rAvg = rSum / dtSum;
			double r2Avg = r2Sum / dtSum;

			rAvgSum += rAvg;
			r2AvgSum += r2Avg;
			rAvg2Sum += rAvg*rAvg;
		}
		double rAvgAvg = rAvgSum / NUMBER_OF_TRIALS;        // <<r>>
		double r2AvgAvg = r2AvgSum / NUMBER_OF_TRIALS;      // <<r^2>>
		double rAvg2Avg = rAvg2Sum / NUMBER_OF_TRIALS;      // <<r>^2>
		double X = r2AvgAvg - rAvgAvg*rAvgAvg;    // <<r^2>> - <<r>>^2
		double Xnew = rAvg2Avg - rAvgAvg*rAvgAvg; // <<r>^2> - <<r>>^2

		// track progress
		std::cout << aRange[a] << " finished\t" << "[" << a+1 << "/" << numPoints << "]\n";
		// write to file
		rvsaFile << std::fixed << std::setprecision(12)
		         << aRange[a] << "\t" << rAvgAvg << "\t" << X << "\t" << Xnew << std::endl;
	}

	return 0;
}
