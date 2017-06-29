#include <iostream>
#include <math.h>

#include "pcg_random.hpp"
#include "topology.hpp"

Topology::Topology(int N, int k, double p)
{
	this->N = N;
	this->k = k;
	this->p = p;

	// connectome is a matrice where rows and columns represent
	// lattice sites. If site i is connected to j then
	// 	connectome[i,j] = connectome[j,i] = 1
	kernelSizes.resize(N, 0);
	regularRing(); // TODO generalize to rewired rings

	// get minimum and maximum number of neighbors
	minNeighbors = maxNeighbors = kernelSizes[0];
	for(int i = 1; i < N; ++i) {
		int size = kernelSizes[i];
		if(size > maxNeighbors) maxNeighbors = size;
		else if(size < minNeighbors) minNeighbors = size;
	}
}

void Topology::printTopology()
{
}

void Topology::regularRing()
{
	// populate the three relevant vectors:
	// kernelSizes - stores the number of neighbors for each site (N elements)
	// kernelList - stores the kernels of each site in sequence (N*k elements)
	// kernelId - stores the begining of each kernel for accessing kernelList (N elements)
	std::vector<bool> connectome(N*N, false);
	int idx;
	for(int i = 0; i < N; ++i) {
		for(int j = 0; j < N; ++j) {
			idx = N*i+j;
			if(distance(i,j) <= k && i != j) { // if i is connected to j on a regular ring
				connectome[idx] = true;
				++kernelSizes[i];
				kernelList.push_back(j);
			}
		}
	}

	// generate kernel indexes
	int sum = 0;
	for(int i = 0; i < N; ++i) {
		kernelId.push_back(sum);
		sum += kernelSizes[i];
	}

	std::cout << "Created regular ring with N=" << N << " and k=" << k << std::endl;
	std::cout << "kernelList vector has " << kernelList.size() << " elements" << std::endl;
}

int Topology::distance(int i, int j)
{
	int distance = abs(i-j);
	if(distance <= N/2) return distance;
	else return N-distance;

	return distance;
}
