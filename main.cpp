#include <iostream>
#include <vector>
#include <random>

#include "pcg_random.hpp"
#include "mylib.hpp"

#define non_deterministic_seed true

int main(int argc, char *argv[]) {

	pcg64 rng(42u, 54u);
	if(non_deterministic_seed) {
		rng.seed(pcg_extras::seed_seq_from<std::random_device>());
	}

	Lattice lattice(10,2);

	return 0;
}
