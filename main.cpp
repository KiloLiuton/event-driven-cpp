#include <iostream>
#include <vector>
#include <random>

#include "pcg_random.hpp"
#include "topology.h"
#include "lattice.h"

#define non_deterministic_seed false

class A {
	int& a;
	pcg64& rng;
public:
	A(int&, pcg64&);

	int roll() {
		return rng(6);
	}
	int getA() { return a; }
	void setA() { a = 666; }
};

A::A(int& pa, pcg64& rng) : a(pa), rng(rng) {}

int main(int argc, char *argv[]) {

	// seed rng
	pcg64 rng(42u, 54u);
	if(non_deterministic_seed) rng.seed(pcg_extras::seed_seq_from<std::random_device>());

	// create a lattice instance
	Lattice lattice(10,3,1.5,0.0,rng);

	std::cout << "\ninitial condition\n";
	lattice.print();

	lattice.step();
	std::cout << "\nafter 1 step\n";
	lattice.print();

	lattice.setCouplingStrength(1.8);
	std::cout << "\nafter setting a=1.8\n";
	lattice.print();

	lattice.step();
	std::cout << "\nafter another step\n";
	lattice.print();

	for(int i = 0; i < 50; ++i) {
		lattice.step();
		lattice.printStates();
	}

	std::cout << "\nafter 50 steps\n";
	lattice.print();

	// get neighbors of site 0
	std::vector<int> foo = lattice.Topology::getNeighbors(0);
	for(const auto& x : foo) std::cout << x << " ";
	std::cout << std::endl;

	lattice.printTopology();

	return 0;
}
