#include <iostream>
#include <math.h>

#include "pcg_random.hpp"
#include "lattice.h"
#include "topology.h"

Lattice::Lattice(int N, int k, double a, double p, pcg64& rng) : Topology(N,k,p), rng(rng) // topology has no default constructor
{
	// set lattice size N, k, and topology at initialization
	this->N = N;
	this->k = k;
	this->a = a;
	this->p = p;

	states.resize(N);
	transitionRates.resize(N);

	initializeStates();
	calculateTransitionsTable();
	updateTransitionRates();
}

void Lattice::initializeStates()
{
	// allocate 'states' vector and randomize its entries
	N0 = N1 = N2 = 0;
	short int state;
	for(int i = 0; i < N; ++i) {
		state = (short int) rng(3);
		states[i] = state;
		switch(state) {
			case 0: ++N0;
					break;
			case 1: ++N1;
					break;
			case 2: ++N2;
					break;
			default: throw std::runtime_error("invalid state value assigned to states vector");
		}
	}
}

int Lattice::expIndex(int k, int dk)
{
	if(k > maxNeighbors || k < minNeighbors || dk > k || dk < -k) {
		throw std::runtime_error("accessing index out of bounds in expTable");
	}
	return (k - minNeighbors) * (minNeighbors + k) + (k + dk);
}

void Lattice::calculateTransitionsTable()
{
	// pre-calculate an exponential table for a particular value of coupling strength 'a'.
	// if there are too many transitions it might be faster to compute the transition as required.
	// transition rate: g = exp[a*(Knext - Ksame)/K]
	// number of possible transitions: (kmax + kmin + 1)*(kmax - kmin + 1)
	transitionsTable.resize((maxNeighbors + minNeighbors + 1) * (maxNeighbors - minNeighbors + 1));
	int i = 0;
	for(int k = minNeighbors; k < maxNeighbors + 1; ++k) {
		for(int ki = -k; ki <= k; ++ki) {
			transitionsTable[i] = exp(a*ki/k);
			++i;
		}
	}
}

void Lattice::updateTransitionRates()
{
	// for every site i get its 'delta' (=knext - ksame) and its number of neighbors 'k'.
	// use this values to access the pre-calculated exp table or calculate it directly.
	for(int i = 0; i < N; ++i) {
		int delta = 0;
		short int currentState = states[i];
		std::vector<int> neighbors = getNeighbors(i);
		for(const auto& s : neighbors) {
			if(states[s] == currentState) --delta;
			else if(states[s] == (currentState+1)%3) ++ delta;
		}
		int k = neighbors.size();

		transitionRates[i] = transitionsTable[expIndex(k,delta)];
	}
}

double Lattice::getOrderParameter()
{
	// calculate the order parameter for the current state
	return sqrt((double) N0*N0 + N1*N1 + N2*N2 - N1*N2 - N0*N1 - N0*N2)/N;
}

void Lattice::print()
{
	std::cout << "states: ";
	for(int i = 0; i < N; ++i) std::cout << states[i] << " ";
	std::cout << std::endl;
	
	std::cout << "populations: " << N0 << " " << N1 << " " << N2 << std::endl;
	std::cout << "min/max neighbors: " << minNeighbors << "," << maxNeighbors << std::endl;

	std::cout << "transition rates: ";
	for(int i = 0; i < N; ++i) std::cout << transitionRates[i] << " ";
	// for(const auto& g : transitionRates) std::cout << g << " ";
	std::cout << std::endl;
}
