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

	void printTopology();

	std::vector<int> kernelList;
	std::vector<int> kernelId;
	std::vector<int> kernelSizes;

private:
	int N, k; // size and number of forward neighbors
	double p; // reconnection probability
	int minNeighbors, maxNeighbors;


	void regularRing();
	int distance(int, int);
};

#endif
