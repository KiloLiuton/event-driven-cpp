#ifndef LATTICE_H_INCLUDED
#define LATTICE_H_INCLUDED

#include <iostream>
#include <vector>
#include <random>

#include "pcg_random.hpp"
#include "topology.h"

class Lattice : public Topology {
public:
	// size, k, coupling strength, rewire probability, pcg64 reference for
	//                                                 a stream of random numbers
	Lattice(int, int, double, double, pcg64&);

	double getOrderParameter();
	void setCouplingStrength(); //TODO
	void print();
private:
	std::vector<short int> states;
	std::vector<int> deltas;
	std::vector<double> transitionRates;
	std::vector<double> transitionsTable;
	int N, k; // size and half the number neighbors
	int N0, N1, N2; // individual state populations
	double p, a; // reconnection probability and coupling strength
	pcg64& rng;
	std::uniform_real_distribution<double> uniform;
	double totalRate;

	void initializeStates();
	void initializeDeltasAndRates();
	int chooseEvent();
	void transitionSite(int);
	int getSiteDelta(int);

	void calculateTransitionsTable();
	int expIndex(int, int);
};

#endif
