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

	Lattice lattice(10,2,1.5,0.0,rng);
	double r = lattice.getOrderParameter();
	lattice.print();
	std::cout << "r = " << r << std::endl;

	std::vector<int> foo = lattice.getNeighbors(0);
	for(const auto& x : foo) std::cout << x << " ";
	std::cout << std::endl;

	lattice.printTopology();

	// rng test
	int a = 6;
	A myClass(a, rng);
	std::cout << "myClass rolls: ";
	pcg64 rng2 = rng; // make a copy of state before rolls
	for(int i = 0; i < 10; ++i) std::cout << myClass.roll() << " ";
	std::cout << std::endl;
	std::cout << "global rolls: ";
	for(int i = 0; i < 10; ++i) std::cout << rng2(6) << " ";
	std::cout << std::endl;

	return 0;
}
