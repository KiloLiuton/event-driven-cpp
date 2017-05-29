#include <iostream>
#include <math.h>

#include "pcg_random.hpp"
#include "lattice.h"
#include "topology.h"

Lattice::Lattice(int N, int k, double couplingStrength, double rewireProb, pcg64& rng) : Topology(N,k,rewireProb), rng(rng), uniform(0.0,1.0)
{
	// set lattice size N, k, and topology at initialization
	this->N = N;
	this->k = k;
	this->couplingStrength = couplingStrength;
	this->rewireProb = rewireProb;
	this->totalRate = 0;

	states.resize(N);
	transitionRates.resize(N);
	deltas.resize(N);
	transitionsTable.resize((Topology::maxNeighbors + Topology::minNeighbors + 1) * (Topology::maxNeighbors - Topology::minNeighbors + 1));

	initializeStates();
	calculateTransitionsTable();
	initializeDeltas();
	initializeRates(); // sets rates and totalRate
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
			case 0:
				++N0;
				break;
			case 1:
				++N1;
				break;
			case 2:
				++N2;
				break;
			default:
					throw std::runtime_error("invalid state value in 'initializeStates'");
		}
	}
}

void Lattice::initializeDeltas()
{
	// get the delta value for each site and get its transition rate
	for(int i = 0; i < N; ++i) {
		std::vector<int> neighbors = Topology::getNeighbors(i);
		int delta = getSiteDelta(i);
		deltas[i] = delta;
	}
}

void Lattice::initializeRates()
{
	totalRate = 0;
	for(int i = 0; i < N; ++i) {
		double g = transitionsTable[expIndex(Topology::kernelSizes[i], deltas[i])];
		transitionRates[i] = g;
		totalRate += g;
	}
}

int Lattice::getSiteDelta(int site)
{
	int delta = 0;
	std::vector<int> neighbors = Topology::getNeighbors(site);
	short int currentState = states[site];
	short int nextState = (currentState+1)%3;
	for(const auto& n : neighbors) {
		short int neighborState = states[n];
		if(neighborState == currentState) --delta;
		else if(neighborState == nextState) ++delta;
	}

	return delta;
}

int Lattice::chooseEvent()
{
	// choose a site to suffer a transition proportionally to its transition rate.
	// an 'event' is the index of the site that suffered a transition
	double partialRate = 0, g = 0;
	double randomRate = uniform(rng) * totalRate;
	for(int event = 0; event < N; ++event) {
		g = transitionRates[event];
		partialRate += g;
		if(randomRate < partialRate) return event;
	}
	throw std::runtime_error("no valid event chosen at function's 'chooseEvent' end");
}

void Lattice::transitionSite(int site)
{
	// this function is called if 'site' transitioned. Then, update its state, delta and transition rate.
	// also updates its neighbors deltas and transition rates.
	short int newState = (states[site]+1)%3;
	short int newNextState = (newState+1)%3;
	states[site] = newState;
	switch(newState) {
		case 0:
			++N0;
			--N2;
			break;
		case 1: 
			++N1;
			--N0;
			break;
		case 2:
			++N2;
			--N1;
			break;
		default:
			throw std::runtime_error("transitioned to an invalid state in 'transitionSite'");
	}

	std::vector<int> neighbors = Topology::getNeighbors(site);
	for(const auto& n : neighbors) {
		short int neighborState = states[n];
		if(neighborState == newState) {
			deltas[site] -= 2;
			deltas[n] -= 1;
		}
		else if(neighborState == newNextState) {
			deltas[site] += 1;
			deltas[n] -= 1;
		}
		else {
			deltas[site] += 1;
			deltas[n] += 2;
		}
		// site 'n' have its delta changed exaclty one time, because only 'site'
		// transitioned and 'n' retains its state
		totalRate -= transitionRates[n];
		double newRate = transitionsTable[expIndex(Topology::kernelSizes[n], deltas[n])];
		transitionRates[n] = newRate;
		totalRate += newRate;
	}
	// 'site' has its delta changed a number of times equal to its kernelSize
	totalRate -= transitionRates[site];
	double newRate = transitionsTable[expIndex(neighbors.size(), deltas[site])];
	transitionRates[site] = newRate;
	totalRate += newRate;
}

void Lattice::calculateTransitionsTable()
{
	// pre-calculate an exponential table for a particular value of coupling strength 'a'.
	// if there are too many transitions it might be faster to compute the transition as required.
	// transition rate: g = exp[a*(Knext - Ksame)/K]
	// number of possible transitions: (kmax + kmin + 1)*(kmax - kmin + 1)
	int i = 0;
	for(int k = minNeighbors; k < maxNeighbors + 1; ++k) {
		for(int ki = -k; ki <= k; ++ki) {
			transitionsTable[i] = exp(couplingStrength*ki/k);
			++i;
		}
	}
}

int Lattice::expIndex(int k, int dk)
{
	// return the one dimensional index with the value of exp(a*dk/k)
	if(k > maxNeighbors || k < minNeighbors || dk > k || dk < -k) {
		throw std::runtime_error("accessing index out of bounds in expTable");
	}
	return (k - minNeighbors) * (minNeighbors + k) + (k + dk);
}

void Lattice::setCouplingStrength(double a)
{
	this->couplingStrength = a;
	calculateTransitionsTable();
	initializeRates();
}

double Lattice::getOrderParameter()
{
	// calculate the order parameter for the current state
	return sqrt((double) N0*N0 + N1*N1 + N2*N2 - N1*N2 - N0*N1 - N0*N2)/N;
}

void Lattice::print()
{
	std::cout << "states: ";
	for(const auto& s : states) std::cout << s << " ";
	std::cout << std::endl;
	
	std::cout << "populations: " << N0 << " " << N1 << " " << N2 << std::endl;
	std::cout << "min/max neighbors: " << minNeighbors << "," << maxNeighbors << std::endl;

	std::cout << "deltas: ";
	for(const auto& d : deltas) std::cout << d << " ";
	std::cout << std::endl;

	std::cout << "transition rates: ";
	for(const auto& g : transitionRates) std::cout << g << " ";
	std::cout << std::endl;
}

