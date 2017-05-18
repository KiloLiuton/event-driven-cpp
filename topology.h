#ifndef TOPOLOGY_H_INCLUDED
#define TOPOLOGY_H_INCLUDED

#include <iostream>
#include <vector>

#include "pcg_random.hpp"

class Topology {
public:
	int minNeighbors, maxNeighbors;

	Topology(int, int, double);
	std::vector<int> getNeighbors(int);
	void printTopology();

private:
	int N, k; // size and number of forward neighbors
	double p; // reconnection probability
	std::vector<bool> connectome;

	void regularRing();
	int index(int, int);
};

#endif
