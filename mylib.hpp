#ifndef MYLIB_INCLUDED
#define MYLIB_INCLUDED

#include <iostream>
#include <vector>
#include <random>

#include "pcg_random.hpp"

class Topology {
public:
	void regularRing(int, int);
	// TODO implement rewiredRing(int, int, double);
	std::vector<int> getNeighbors(int);

	void print();
private:
	std::vector<bool> connectome;
	int index(int, int);
	int N, k;
};

class Lattice {
public:
	Lattice(int, int);
	// TODO implement rewiring: Lattice(int, int, double);
	Topology topology;

	double getOrderParameter();
	void initializeStates(pcg64& rng);
	void updateTransitionRates();

	void print();
private:
	std::vector<short int> states;
	std::vector<double> transitionRates;
	int N;
	int N0 = 0, N1 = 0, N2 = 0;
	int k;
};

#endif
