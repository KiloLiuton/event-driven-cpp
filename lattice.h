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
	int getPop(short int);
	void step();
	void setCouplingStrength(double);
	void print();
	void printStates();

private:
	std::vector<short int> states;
	std::vector<int> deltas;
	std::vector<double> transitionRates, transitionsTable;
	double totalRate, rewireProb, couplingStrength;
	int N, N0, N1, N2, k; // size, populations and half the number neighbors
	pcg64& rng;
	std::uniform_real_distribution<double> uniform;

	void initializeStates();
	void initializeDeltas();
	void initializeRates();
	void calculateTransitionsTable();
	void transitionSite(int);
	int chooseEvent();
	int getSiteDelta(int);
	int expIndex(int, int);
};

#endif
