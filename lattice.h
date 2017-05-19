#ifndef LATTICE_H_INCLUDED
#define LATTICE_H_INCLUDED

#include <iostream>
#include <vector>

#include "pcg_random.hpp"
#include "topology.h"

class Lattice : public Topology {
public:
	Lattice(int, int, double, pcg64&);
	// TODO implement rewiring: Lattice(int, int, double);

	void initializeStates(pcg64& rng);
	void updateTransitionRates(); // TODO
	double getOrderParameter();

	void print();
private:
	std::vector<short int> states;
	std::vector<double> transitionRates;
	std::vector<double> transitionsTable;

	void calculateTransitionsTable(); // TODO
	int N, k; // size and number of forward neighbors
	int N0, N1, N2; // individual state populations
	double p, a; // reconnection probability and coupling strength
};

#endif
