#include <iostream>
#include <math.h>

#include "pcg_random.hpp"
#include "topology.hpp"

Topology::Topology(int N, int k, double p)
{
	this->N = N;
	this->k = k;
	this->p = p;

	connectome.resize((N*N-N)/2, false);
	kernelSizes.resize(N, 0);
	regularRing(); // TODO generalize to rewired rings

	// get minimum and maximum number of neighbors
	int size = getNeighbors(0).size();
	minNeighbors = maxNeighbors = size;
	for(int i = 1; i < N; ++i) {
		size = getNeighbors(i).size();
		if(size > maxNeighbors) maxNeighbors = size;
		else if(size < minNeighbors) minNeighbors = size;
	}
}

std::vector<int> Topology::getNeighbors(int x)
{
	// return all neighbors of site x by sweeping row and then column
	std::vector<int> neighbors;

	for(int j = x+1; j < N; ++j) {
		if(connectome[index(x,j)]) {
			neighbors.push_back(j);
		}
	}
	for(int i = 0; i < x; ++i) {
		if(connectome[index(i,x)]) {
			neighbors.push_back(i);
		}
	}

	return neighbors;
}

void Topology::printTopology()
{
	// print the half-matrix connectome
	for(int i = 0; i < N-1; ++i) {
		for(int j = N-i; j < N; ++j) { std::cout << "- "; }
		for(int j = i+1; j < N; ++j) { std::cout << connectome[index(i,j)] << " "; }
		std::cout << std::endl;
	}
}

void Topology::regularRing()
{
	// allocate connectome with a regular ring topology
	int idx;
	for(int i = 0; i < N-1; ++i) {
		for(int j = i+1; j < N; ++j) {
			idx = index(i,j);
			if(j <= i+k || j >= N-k+i) { // if i is connected to j increment its kernels
				connectome[idx] = true;
				++kernelSizes[i];
				++kernelSizes[j];
			}
		}
	}
	std::cout << "Created regular ring with N=" << N << " and k=" << k << std::endl;
}



int Topology::index(int i, int j)
{
	// convert (i,j) notation to compact notation with bounds checking
	if(i >= N-1 || i < 0 || j < 1 || j >= N) {
		throw std::runtime_error("index out of bounds in topology");
	}
	int idx;
	idx = i*(N-2) - (i*i-i)/2 + j - 1;
	return idx;
}
