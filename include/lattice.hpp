#ifndef LATTICE_H_INCLUDED
#define LATTICE_H_INCLUDED

#include <iostream>
#include <vector>
#include <random>

#include "pcg_random.hpp"
#include "topology.hpp"

class Lattice : public Topology {
public:
	// size, k, coupling strength, rewire probability, pcg64 reference for
	//                                                 a stream of random numbers
	Lattice(int const, int const, double const, double, pcg64&);

	double getOrderParameter();
	int getPop(short int);
	double step();
	void reset();
	void resetToCoupling(double);
	void setCouplingStrength(double);
	void print();
	void printStates();
	void printPops();
	size_t relaxationRun(int trail, double threshold, const size_t MAX_ITERS, std::ofstream& outputFile);

private:
	const int N; // size and neighbors
	std::vector<short int> states;
	std::vector<int> deltas;
	std::vector<double> transitionRates, transitionsTable;
	double totalRate, couplingStrength;
	int N0, N1, N2; // populations
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
