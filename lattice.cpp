#include <iostream>
#include <math.h>

#include "pcg_random.hpp"
#include "lattice.h"
#include "topology.h"

// set lattice size N, k, and topology at initialization
Lattice::Lattice(int N, int k, double p, pcg64& rng) : Topology(N,k,p) // topology has no default constructor
{
	this->N = N;
	this->k = k;
	this->p = p;
	states.resize(N);
	transitionRates.resize(N);

	initializeStates(rng);
}

// allocate 'states' vector and randomize its entries
void Lattice::initializeStates(pcg64& rng)
{
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

void Lattice::calculateTransitionsTable()
{
	// g = exp[a*(Knext - Ksame)/K]
}

void Lattice::updateTransitionRates()
{

}

// calculate the order parameter for the current state
double Lattice::getOrderParameter()
{
	return sqrt((double) N0*N0 + N1*N1 + N2*N2 - N1*N2 - N0*N1 - N0*N2)/N;
}

void Lattice::print()
{
	std::cout << "states: ";
	for(int i = 0; i < N; ++i) { std::cout << states[i] << " ";	}
	std::cout << std::endl;
	std::cout << "populations: ";
	std::cout << N0 << " " << N1 << " " << N2;
	std::cout << std::endl;
}
