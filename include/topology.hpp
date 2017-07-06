#ifndef TOPOLOGY_H_INCLUDED
#define TOPOLOGY_H_INCLUDED

#include <iostream>
#include <vector>

#include "pcg_random.hpp"

class Topology {
public:
	Topology(const int, const int, const double); // constructor

	int getMaxNeighbors() const { return maxNeighbors; }
	int getMinNeighbors() const { return minNeighbors; }

	void printTopology() const; // graphically print connectivity matrix
	void printKernels() const; // print kernels as lists of indexes

	std::vector<int> kernelList; // store all kernels sequentially
	std::vector<int> kernelId; // store an index to the begining of each kernel
	std::vector<int> kernelSizes; // store all kernel sizes

private:
	const int N, k; // size and number of forward neighbors
	const double p; // reconnection probability
	int minNeighbors, maxNeighbors;

	void createRing();
	void printKernel(int) const;
	bool isInKernel(int, int) const;
};

#endif
