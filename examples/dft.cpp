#include "../efgen/dft.hpp"
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

// Example DFT usage: multiplying polynomials.

using complex = std::complex<double>;

int main() {
	int n;
	std::cin >> n;
	
	std::vector<int> coeffs_a(n);
	for (int i = 0; i < n; ++i) {
		std::cin >> coeffs_a[i];
	}
	
	int m;
	std::cin >> m;
	
	std::vector<int> coeffs_b(m);
	for (int i = 0; i < m; ++i) {
		std::cin >> coeffs_b[i];
	}
	
	int N = 2 << (std::numeric_limits<unsigned>::digits - __builtin_clz(std::max(n, m) - 1));
	
	std::vector<complex> A(N);
	std::vector<complex> B(N);
	for (int i = 0; i < n; ++i) A[i] = coeffs_a[i];
	for (int i = 0; i < m; ++i) B[i] = coeffs_b[i];
	
	complex omega = std::polar(1.0, 2.0 * M_PI / N);
	
	dft(N, A.data(), omega);
	dft(N, B.data(), omega);
	
	for (int i = 0; i < N; ++i) {
		A[i] *= B[i];
	}
	
	dft(N, A.data(), 1.0 / omega);
	
	std::vector<int> coeffs_ab(N);
	for (int i = 0; i < N; ++i) {
		coeffs_ab[i] = lround(A[i].real() / N);
	}
	
	while (coeffs_ab.size() > 1 && coeffs_ab.back() == 0) coeffs_ab.pop_back();
	
	std::cout << "Product:";
	for (int coeff : coeffs_ab) std::cout << ' ' << coeff;
	std::cout << '\n';
}
