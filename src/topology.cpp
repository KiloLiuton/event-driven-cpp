#include <iostream>
#include <math.h>
#include <random>

#include "pcg_random.hpp"
#include "topology.hpp"

Topology::Topology(
		int const N,
		int const k,
	   	double const p,
		bool const NON_DETERMINISTIC_TOPOLOGY
		) : N(N), k(k), p(p), NON_DETERMINISTIC_TOPOLOGY(NON_DETERMINISTIC_TOPOLOGY)
{
	createRing();

	// get minimum and maximum number of kernel sizes
	minNeighbors = maxNeighbors = kernelSizes[0];
	for (int i = 1; i < N; ++i) {
		int size = kernelSizes[i];
		if(size > maxNeighbors) maxNeighbors = size;
		else if(size < minNeighbors) minNeighbors = size;
	}
}

void Topology::createRing()
{
	// populate the three relevant vectors:
	// kernelSizes - store the number of edges emanating from each vertex (N elements)
	// kernelList - stores a list of connected vertices (kernel) for each vertex (N*k elements)
	// kernelId - stores the begining of kernel for each vertex (N elements)
	//
	// Start by making a regular ring with trivial kernel sizes of 2k for every vertex
	for (int i = 0; i < N; ++i) {
		for (int j = -k; j <= k; ++j) {
			if (j == 0) continue; // a site cannot have itself in its own kernel
			int idx = i + j;
			if (idx < 0) idx += N;
			else if (idx >= N) idx -= N;
			kernelList.push_back(idx);
		}
		kernelSizes.push_back(2*k);
	}

	// generate kernel indexes
	int sum = 0;
	for(int i = 0; i < N; ++i) {
		kernelId.push_back(sum);
		sum += kernelSizes[i];
	}

	if (p == 0.0) {
		std::cout << "\nCreated regular ring with N=" << N << " and k=" << k << "\n";
	} else {
		// make an rng
		pcg64 rng(42u, 54u);
		std::uniform_real_distribution<double> uniform(0.0, 1.0);
		if(NON_DETERMINISTIC_TOPOLOGY) rng.seed(pcg_extras::seed_seq_from<std::random_device>());

		// for each vertex, loop only through the k clockwise edges
		for (int currentVertex = 0; currentVertex < N-1; ++currentVertex) {
			for (int j = 1; j <= k; j++) {
				if (uniform(rng) > p) continue; // rewire edge with probability p

				int cutVertex = currentVertex + j;
				if (cutVertex >= N) cutVertex -= N;

				// prevent rewiring from leaving isolated vertices and also chosing invalid edges
				if (kernelSizes[cutVertex] <= 1) continue;
				int randomVertex = rng(N);
				while (isInKernel(currentVertex, randomVertex) || randomVertex == currentVertex) randomVertex = rng(N);

				// swap the endpoint of the edge
				std::vector<int>::iterator first = kernelList.begin() + kernelId[currentVertex];
				std::vector<int>::iterator last = first + kernelSizes[currentVertex];
				std::vector<int>::iterator swapPoint = std::find(first, last, cutVertex);
				if (swapPoint != last) *swapPoint = randomVertex;
				else throw std::runtime_error("cut vertex not in current vertex kernel");

				// remove the current vertex from the cut vertex kernel
				first = kernelList.begin() + kernelId[cutVertex];
				last = first + kernelSizes[cutVertex];
				std::vector<int>::iterator removePoint = std::find(first, last, currentVertex);
				if (removePoint != last) kernelList.erase(removePoint);
				else throw std::runtime_error("Trying to remove inexistent vertex from kernel");
				for (int i = cutVertex + 1; i < N; ++i) kernelId[i]--;
				kernelSizes[cutVertex]--;

				// insert the current vertex at the end of randomly selected vertex kernel
				std::vector<int>::iterator insertPoint = kernelList.begin() + kernelId[randomVertex] + kernelSizes[randomVertex];
				kernelList.insert(insertPoint, currentVertex);
				// increment indexes after the insertion point
				for (int i = randomVertex + 1; i < N; ++i) kernelId[i]++;
				kernelSizes[randomVertex]++;
			}
		}
		std::cout << "\nCreated rewired ring with N="<<N<<" k="<<k<<" and p="<<p<<"\n";
	}
}

bool Topology::isInKernel(int kernelNum, int element) const
{
	// return if element is in kernel number kernelNum
	std::vector<int>::const_iterator start = kernelList.begin() + kernelId[kernelNum];
	std::vector<int>::const_iterator finish = start + kernelSizes[kernelNum];

	if (std::find(start, finish, element) != finish) return true;
	return false;
}

void Topology::printKernel(int i) const
{
	int first = kernelId[i];
	int last = first + kernelSizes[i];
	for (int i = first; i < last; ++i) std::cout << kernelList[i] << " ";
	std::cout << "\n";
}

void Topology::printTopology() const
{
	// print a square connectivity matrix
	std::cout << "\\ ";
	for(int i = 0; i < N; ++i) std::cout << i << " ";
	std::cout << "\n";
	for(int i = 0; i < N; ++i) {
		std::cout << i << " ";
		std::vector<int>::const_iterator start = kernelList.begin() + kernelId[i];
		std::vector<int>::const_iterator end = kernelList.begin() + kernelId[i] + kernelSizes[i];
		for (int j = 0; j < N; ++j) {
			if (std::find(start, end, j) != end) {
				std::cout << "* ";
			} else if (i == j) {
				std::cout << "\\ ";
			} else {
				std::cout << ". ";
			}
		}
		std::cout << "\n";
	}
}

void Topology::printKernels() const
{
	for (int i = 0; i < N; ++i) {
		printKernel(i);
		std::cout << " ";
	}
	std::cout << "\n";
}

void Topology::printToFile() const
{
	std::cout << "print to file not implemented yet\n";
}
