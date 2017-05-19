#include <iostream>
#include <vector>
#include <random>

#include "pcg_random.hpp"
#include "topology.h"
#include "lattice.h"

#define non_deterministic_seed true

int main(int argc, char *argv[]) {

	// seed rng
	pcg64 rng(42u, 54u);
	if(non_deterministic_seed) rng.seed(pcg_extras::seed_seq_from<std::random_device>());

	Lattice lattice(10,2,0.0,rng);
	double r = lattice.getOrderParameter();
	lattice.print();
	std::cout << "r = " << r << std::endl;

	std::vector<int> foo = lattice.getNeighbors(0);
	for(const auto& x : foo) std::cout << x << " ";
	std::cout << std::endl;

	lattice.printTopology();

	return 0;
}
