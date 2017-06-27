#ifndef TOPOLOGY_H_INCLUDED
#define TOPOLOGY_H_INCLUDED

#include <iostream>
#include <vector>

#include "pcg_random.hpp"

class Topology {
public:
	Topology(int, int, double); // constructor

	int getMaxNeighbors() { return maxNeighbors; }
	int getMinNeighbors() { return minNeighbors; }
	int getKernelSize(int i) { return kernelSizes[i]; }
	std::vector<int> getNeighbors(int);

	void printTopology();

private:
	int N, k; // size and number of forward neighbors
	double p; // reconnection probability
	int minNeighbors, maxNeighbors;

	std::vector<bool> connectome;
	std::vector<int> kernelSizes;

	void regularRing();
	int index(int, int);
	int distance(int, int);
};

#endif
