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

	// seed rng
	pcg64 rng(42u, 54u);
	if(non_deterministic_seed) rng.seed(pcg_extras::seed_seq_from<std::random_device>());

	// define lattice parameters:
	const int SIZE = 12;
	const int K = 2;
	const double REWIRE_PROBABILITY = 0.0; // TODO : this is not implemented yet! only regular rings are created
	double couplingStrength = 1.6;

	// create a lattice instance
	Lattice lattice(SIZE, K, REWIRE_PROBABILITY, couplingStrength, rng);

	// print topology and initial condition FIXME: (may print incorrectly if SIZE is too big)
	lattice.printTopology();
	std::cout << "\nINITIAL CONDITION\n";
	lattice.print();

	// set simulation parameters
	const int ITERS = 2000;
	const int TRIALS = 1;

	// prepare output file
	std::ostringstream oss;
	std::string dataPath ("data/");
	oss << "N=" << SIZE << "k=" << K << "p=" << REWIRE_PROBABILITY << "ITERS=" << ITERS << "TRIALS=" << TRIALS << ".txt";
	std::string filename = oss.str();

	std::ofstream orderParameterFile (dataPath + filename);
	if(!orderParameterFile.is_open()) throw std::runtime_error("failed to open file");

	// run system dynamics
	orderParameterFile << lattice.getOrderParameter() << std::endl;
	lattice.step();
	orderParameterFile << lattice.getOrderParameter() << std::endl;

	return 0;
}
