#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <algorithm>


#include "pcg_random.hpp"
#include "lattice.hpp"
#include "topology.hpp"

Lattice::Lattice(int N, int k, double p, double couplingStrength, pcg64& rng) : Topology(N,k,p), N(N), rng(rng), uniform(0.0,1.0)
{
	// set lattice size N, k, and topology at initialization
	this->couplingStrength = couplingStrength;
	this->totalRate = 0;

	states.resize(N);
	transitionRates.resize(N);
	deltas.resize(N);
	int max = Topology::getMaxNeighbors();
	int min = Topology::getMinNeighbors();
	// resize transitions table to accomodate all possible transition values
	transitionsTable.resize((max + min + 1) * (max - min + 1));

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
		deltas[i] = getSiteDelta(i);
	}
}

void Lattice::initializeRates()
{
	totalRate = 0;
	for(int i = 0; i < N; ++i) {
		double g = transitionsTable[expIndex(Topology::getKernelSize(i), deltas[i])];
		transitionRates[i] = g;
		totalRate += g;
	}
}

int Lattice::getSiteDelta(int site)
{
	int delta = 0;
	short int currentState = states[site];
	short int nextState = (currentState+1)%3;
	std::vector<int> neighbors = Topology::getNeighbors(site);
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
	throw std::runtime_error("no valid event chosen at function 'chooseEvent's end");
}

void Lattice::transitionSite(int site)
{
	// this function is called if 'site' transitioned. Then, update its state, delta and transition rate.
	// also updates its neighbors deltas and transition rates.

	// update site state and populations
	short int newState = (states[site]+1)%3;
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

	// update neighbors states and all deltas
	//	'site' has its delta changed a number of times equal to its kernelSize
	//	each 'n' retains its state and have its delta changed exaclty one time
	short int newNextState = (newState+1)%3;
	for(const auto& n : Topology::getNeighbors(site)) {
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
		double newRate = transitionsTable[expIndex(Topology::getKernelSize(n), deltas[n])];
		totalRate += newRate;
		totalRate -= transitionRates[n];
		transitionRates[n] = newRate;
	}
	double newRate = transitionsTable[expIndex(Topology::getKernelSize(site), deltas[site])];
	totalRate += newRate;
	totalRate -= transitionRates[site];
	transitionRates[site] = newRate;
}

void Lattice::calculateTransitionsTable()
{
	// pre-calculate an exponential table for a particular value of coupling strength 'a'.
	// if there are too many transitions it might be faster to compute the transition as required.
	// transition rate: g = exp[a*(Knext - Ksame)/K]
	// number of possible transitions: (kmax + kmin + 1)*(kmax - kmin + 1)
	int i = 0;
	for(int k = Topology::getMinNeighbors(); k < Topology::getMaxNeighbors() + 1; ++k) {
		for(int ki = -k; ki <= k; ++ki) {
			transitionsTable[i] = exp(couplingStrength*ki/k);
			++i;
		}
	}
}

int Lattice::expIndex(int k, int dk)
{
	// use this function to get the correct value of the exponential for k and dk.
	// return the one dimensional index with the value of exp(a*dk/k).
	if(k > Topology::getMaxNeighbors() || k < Topology::getMinNeighbors() || dk > k || dk < -k) {
		throw std::runtime_error("accessing index out of bounds in expTable");
	}
	return (k - Topology::getMinNeighbors()) * (Topology::getMinNeighbors() + k) + (k + dk);
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

double Lattice::step()
{
	// this function runs the model dynamics for one step. In the evet driven paradigm this
	// means that one event will occur for every call of this function, regardless of the time
	// elapsed.
	// return value is the expected time this state will last until next transition.
	int event = chooseEvent();
	transitionSite(event);
	double expectedTime = 1/totalRate;

	return expectedTime;
}

void Lattice::reset()
{
	initializeStates();
	initializeDeltas();
	initializeRates();
}

void Lattice::resetToCoupling(double a)
{
	initializeStates();
	initializeDeltas();
	setCouplingStrength(a);
}

int Lattice::getPop(short int state)
{
	switch(state) {
		case 0:
			return N0;
		case 1:
			return N1;
		case 2:
			return N2;
		default:
			throw std::runtime_error("invalid state queried at 'getPop'");
	}
}

size_t Lattice::relaxationRun(int trail, double threshold, const size_t MAX_ITERS, std::ofstream& file)
{
	// run a single trial to determine relaxation. Relaxation is found when
	//    the average order parameter doesn't change more than threshold
	//    for trail steps.

	// start by running 'trail' steps and storing the 'r' values.
	double highestAvg = 0; // minimum possible value of r is 0
	double lowestAvg = 1; // maximum possible value of r is 1
	double totalTime = 0;
	std::vector<double> trackR(trail);
	size_t relaxationPeriod = 1;
	trackR[0] = getOrderParameter();
	for(int i = 1; i < trail; ++i) {
		++relaxationPeriod;
		double dt = step();
		double r = getOrderParameter();
		trackR[i] = r;
		totalTime += dt;

		file << std::fixed << std::setprecision(12)
		     << totalTime << "\t" << r << std::endl;
	}
	// get average of the first block of 'trail' events
	double avg = std::accumulate(trackR.begin(), trackR.end(), 0) / trackR.size();
	highestAvg = std::max(highestAvg, avg);
	lowestAvg = std::min(lowestAvg, avg);

	// for each next step, check if the new average is in an interval of 'threshold'
	//     around the previous average. If this happens 'trail' consecutive times, break.
	int count = 0;
	for(size_t i = trail; i < MAX_ITERS; ++i) {
		++relaxationPeriod;
		double dt = step();
		double r = getOrderParameter();
		double nextAvg = avg + (r - trackR[0]) / trail;
		trackR.erase(trackR.begin());
		trackR.push_back(r);
		totalTime += dt;
		if(std::fabs(avg - nextAvg) < threshold) ++count;
		else count = 0;
		avg = nextAvg;
		highestAvg = std::max(highestAvg, avg);
		lowestAvg = std::min(lowestAvg, avg);

		file << std::fixed << std::setprecision(12)
		     << totalTime << "\t" << r << "\t" << getPop(0) << "\t" << getPop(1) << std::endl;

		if(count > trail) break;
	}

	return relaxationPeriod;
}

void Lattice::print()
{
	std::cout << "states: ";
	for(const auto& s : states) std::cout << s << " ";
	std::cout << std::endl;

	std::cout << "deltas: ";
	for(const auto& d : deltas) std::cout << d << " ";
	std::cout << std::endl;
	
	std::cout << "populations: " << N0 << " " << N1 << " " << N2 << std::endl;
	std::cout << "min/max neighbors: " << Topology::getMinNeighbors() << "," << Topology::getMaxNeighbors() << std::endl;

	std::cout << "transition rates: ";
	std::cout.precision(3);
	for(const auto& g : transitionRates) std::cout << g << " ";
	std::cout << std::endl;

	std::cout << "r = " << getOrderParameter() << std::endl;
}

void Lattice::printStates()
{
	for(auto s : states) std::cout << s << " ";
}

void Lattice::printPops()
{
	std::cout << N0 << " " << N1 << " " << N2;
}
