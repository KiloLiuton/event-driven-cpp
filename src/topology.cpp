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

int Topology::index(int i, int j) { return N*i+j; }

std::vector<int> Topology::getNeighbors(int x)
{
	// return all neighbors of site x by sweeping row
	std::vector<int> neighbors;

	for(int j = 0; j < N; ++j) {
		if(connectome[index(x,j)]) {
			neighbors.push_back(j);
		}
	}

	return neighbors;
}

void Topology::printTopology()
{
	// print the connectome
	for(int i = 0; i < N; ++i) {
		for(int j = 0; j < N; j++) {
			std::cout << connectome[index(i,j)] << " ";
		}
		std::cout << std::endl;
	}
}

void Topology::regularRing()
{
	// populate connectome with regular ring connectivity
	connectome.resize(N*N, false);
	int idx;
	for(int i = 0; i < N; ++i) {
		for(int j = 0; j < N; ++j) {
			idx = index(i,j);
			if(distance(i,j) <= k && i != j) { // if i is connected to j increment its kernels
				connectome[idx] = true;
				++kernelSizes[i];
				++kernelSizes[j];
			}
		}
	}
	std::cout << "Created regular ring with N=" << N << " and k=" << k << std::endl;
}

int Topology::distance(int i, int j)
{
	int distance = abs(i-j);
	if(distance <= N/2) return distance;
	else return N-distance;

	return distance;
}
