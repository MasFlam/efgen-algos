//*==[efgen-dft]==========================*//
// Fairly efficient but still generic DFT. //
// - author: Łukasz Drukała                //
// - version: 2024-12-07                   //
// Released to the public domain.          //
// Attribution would be welcome though ;)  //
//*=======================================*//
#pragma once

#include <cassert>
#include <cstddef>
#include <limits>
#include <utility>

template<typename T>
constexpr void dft_permute(size_t k, T A[]) noexcept {
	size_t n = size_t(1) << k;
	for (size_t i = 0; i < n; ++i) {
		// Reverse the k least significant bits.
		size_t j = 0;
		for (int b = 0; b < k; ++b) {
			j = (j << 1) | ((i >> b) & 1);
		}
		// Only swap once.
		if (i < j) std::swap(A[i], A[j]);
	}
}

// To calculate the inverse DFT use dft(n, A, omega^-1) and divide each element by n.
template<typename T>
constexpr void dft(size_t n, T A[], T omega) noexcept {
	if (n <= 1) return;
	
	// k = log2(n)
	size_t k = std::numeric_limits<unsigned long long>::digits - __builtin_clzll(n - 1);
	
	// Assert power-of-two length.
	assert(n == size_t(1) << k);
	
	// Permute the data array beforehand. This gives us:
	// - better cache efficiency in the main algorithm,
	// - in-place computation.
	dft_permute(k, A);
	
	// Precompute omega^2^(k-lvl) for each lvl = 1, ..., k.
	T omegas[k+1];
	omegas[k] = omega;
	for (size_t lvl = k-1; lvl > 0; --lvl) {
		omegas[lvl] = omegas[lvl+1] * omegas[lvl+1];
	}
	
	// Simulate the recursion.
	for (size_t lvl = 1; lvl <= k; ++lvl) {
		size_t blocksz = size_t(1) << lvl;
		for (size_t i = 0; i < n; i += blocksz) {
			T ompow = 1;
			for (size_t j = 0; j < blocksz/2; ++j) {
				T a = A[i + j];
				T b = ompow * A[i + blocksz/2 + j];
				A[i + j]             = a + b;
				A[i + blocksz/2 + j] = a - b;
				ompow *= omegas[lvl];
			}
		}
	}
}
