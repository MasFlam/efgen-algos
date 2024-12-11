#include "../efgen/dft.hpp"
#include <cassert>
#include <iostream>
#include <vector>
using ull = unsigned long long;

// Example DFT usage: multiplying polynomials. (DFT over a finite field)

// MOD^2 would not fit in a signed int64, so let's use unsigned since we wrap these in struct modint anyway.
constexpr ull MOD = (1ull << 30) * 3 + 1; // 2 and 3 are the only divisors of phi(MOD), so finding a generator is easy
constexpr ull GEN = 5;                    // ord(g) = MOD-1
constexpr ull OMEGA = 125;                // omega = g^3 = primitive root of unity of order 2^30

// Returns the primitive root of unity of order 2^logn.
constexpr ull OMEGA_FOR(int logn) noexcept {
	assert(0 <= logn && logn <= 30);
	ull omega = OMEGA;
	for (int i = logn; i < 30; ++i) {
		omega = omega * omega % MOD;
	}
	return omega;
}

constexpr ull modpow(ull x, ull n) noexcept {
	if (n == 0) return 1;
	if (n == 1) return x;
	if (n % 2 == 0) {
		ull r = modpow(x, n/2);
		return r * r % MOD;
	} else {
		ull r = modpow(x, n-1);
		return r * x % MOD;
	}
}

struct modint {
	ull val{};
	constexpr modint() noexcept = default;
	constexpr modint(ull x) noexcept : val(x % MOD) {};
	constexpr operator ull() const noexcept { return val; }
};
// Calculates the inverse using Fermat's little theorem.
constexpr modint inv(modint x) noexcept { return {modpow(x, MOD-2)}; }
constexpr modint operator+(modint x, modint y) noexcept { return {(x.val + y.val) % MOD}; }
constexpr modint operator-(modint x, modint y) noexcept { return {(x.val - y.val + MOD) % MOD}; }
constexpr modint operator*(modint x, modint y) noexcept { return {(x.val * y.val) % MOD}; }
constexpr modint operator/(modint x, modint y) noexcept { return {(x.val * modpow(y.val, MOD-2)) % MOD}; }
constexpr modint operator+(modint x) noexcept { return x; }
constexpr modint operator-(modint x) noexcept { return {MOD - x.val}; }
constexpr modint &operator+=(modint &x, modint y) noexcept { return x = x + y; }
constexpr modint &operator-=(modint &x, modint y) noexcept { return x = x - y; }
constexpr modint &operator*=(modint &x, modint y) noexcept { return x = x * y; }
constexpr modint &operator/=(modint &x, modint y) noexcept { return x = x / y; }
std::ostream &operator<<(std::ostream &stream, modint x) { return stream << x.val; }

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
	
	int logN = 1 + (std::numeric_limits<unsigned>::digits - __builtin_clz(std::max(n, m) - 1));
	int N = 1 << logN;
	
	std::vector<modint> A(N);
	std::vector<modint> B(N);
	for (int i = 0; i < n; ++i) A[i] = coeffs_a[i];
	for (int i = 0; i < m; ++i) B[i] = coeffs_b[i];
	
	modint omega = OMEGA_FOR(logN);
	
	for (int i = 0; i < N; ++i) {
		std::cout << "A[i], B[i] = " << A[i] << ", " << B[i] << '\n';
	}
	std::cout << "omega = " << omega << '\n';
	
	dft(N, A.data(), omega);
	dft(N, B.data(), omega);
	
	for (int i = 0; i < N; ++i) {
		std::cout << "A[i], B[i] = " << A[i] << ", " << B[i] << '\n';
		A[i] *= B[i];
	}
	
	dft(N, A.data(), inv(omega));
	
	modint invN = inv(modint(N));
	std::vector<int> coeffs_ab(N);
	for (int i = 0; i < N; ++i) {
		coeffs_ab[i] = A[i] * invN;
	}
	
	while (coeffs_ab.size() > 1 && coeffs_ab.back() == 0) coeffs_ab.pop_back();
	
	std::cout << "Product:";
	for (int coeff : coeffs_ab) std::cout << ' ' << coeff;
	std::cout << '\n';
}
