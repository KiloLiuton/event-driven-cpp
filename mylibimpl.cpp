#include <iostream>
#include <random>
#include <string>

#include "pcg_random.hpp"
#include "mylib.hpp"

// allocate connectome with a regular ring topology
void Topology::regularRing(int size, int connectivity)
{
	N = size;
	k = connectivity;
	connectome.resize((N*N-N)/2, false);
	int idx;
	for(int i = 0; i < N-1; ++i) {
		for(int j = i+1; j < N; ++j) {
			idx = index(i,j);
			if(j <= i+k || j >= N-k+i) {
				connectome[idx] = true;
			}
		}
	}
	std::cout << "Created regular ring with N=" << N << " and k=" << k << std::endl;
}

// return all neighbors of site x by sweeping row and then column
std::vector<int> Topology::getNeighbors(int x)
{
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

// print the half-matrix connectome
void Topology::print()
{
	for(int i = 0; i < N-1; ++i) {
		for(int j = N-i; j < N; ++j) std::cout << "- ";
		for(int j = i+1; j < N; ++j) {
			std::cout << connectome[index(i,j)] << " ";
		}
		std::cout << std::endl;
	}
}

// convert (i,j) notation to compact notation with bounds checking
int Topology::index(int i, int j)
{
	if(i >= N-1 || i < 0 || j < 1 || j >= N) {
		throw std::runtime_error("accessing index out of bounds");
	}
	int idx;
	idx = i*(N-2) - (i*i-i)/2 + j - 1;
	return idx;
}

// initialize the lattice with N sites and 2*k neighbors per site
Lattice::Lattice(int size, int neighbors)
{
	N = size;
	k = neighbors;
	topology.regularRing(N, k);
}

// allocate states vector and randomize its entries
void Lattice::initializeStates(pcg64& rng)
{
	states.resize(N);
	int state;
	for(int i = 0; i < N; ++i) {
		state = rng(3);
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
