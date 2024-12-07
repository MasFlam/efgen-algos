#include "../efgen/segtree.hpp"
using ll = long long;

// Example traits: a tree with "add x on segment" and "get sum of segment" operations.
struct MyTraits {
	using Val = std::pair<ll, int>;
	using Mod = ll;
	static Val neutral() { return {0, 1}; }
	static Mod ident() { return 0; }
	static Val join(const Val &a, const Val &b) { return {a.first + b.first, a.second + b.second}; }
	static Mod compose(const Mod &f, const Mod &g) { return f + g; }
	static Val apply(const Mod &f, const Val &x) { return {x.first + f * x.second, x.second}; }
};

int main() {
	std::pair<ll, int> data[] {
		{2, 1},
		{1, 1},
		{3, 1},
		{7, 1},
	};
	SegTree<MyTraits> tree(from_iter_tag{}, std::begin(data), std::end(data));
}
