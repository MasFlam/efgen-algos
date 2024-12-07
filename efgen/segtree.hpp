//*==[efgen-segtree]===========================*//
// Efficient and generic segment tree template. //
// - author: Łukasz Drukała                     //
// - version: 2024-12-06                        //
// Released to the public domain.               //
// Attribution would be welcome though ;)       //
//*============================================*//
#pragma once

#include <iterator>
#include <limits>
#include <vector>

// See SegTree's constructors.
struct from_iter_tag {};

template<typename TraitsT>
class SegTree {
public:
	using Val = typename TraitsT::Val;
	using Mod = typename TraitsT::Mod;
private:
	static Val neutral() { return TraitsT::neutral(); }
	static Mod ident() { return TraitsT::ident(); }
	static Val join(const Val &x, const Val &y) { return TraitsT::join(x, y); }
	static Mod compose(const Mod &f, const Mod &g) { return TraitsT::compose(f, g); }
	static Val apply(const Mod &f, const Val &x) { return TraitsT::apply(f, x); }
	
	struct Node {
		Val val{neutral()};
		Mod mod{ident()};
		
		Val modval() const {
			return apply(mod, val);
		}
	};
	
	constexpr static size_t calc_nleaves(size_t n) noexcept {
		if (n <= 4) return 4;
		return size_t(1) << (std::numeric_limits<unsigned long long>::digits - __builtin_clzll(n - 1));
	}
	
	size_t m_length;
	size_t m_nleaves;
	std::vector<Node> m_data;
	
	void recalc(size_t node) {
		m_data[node].val = join(m_data[2 * node + 0].modval(), m_data[2 * node + 1].modval());
	}
	
	void downprop(size_t node, bool update_val) {
		if (update_val) m_data[node].val = m_data[node].modval();
		m_data[2 * node + 0].mod = compose(m_data[node].mod, m_data[2 * node + 0].mod);
		m_data[2 * node + 1].mod = compose(m_data[node].mod, m_data[2 * node + 1].mod);
		m_data[node].mod = ident();
	}
	
public:
	
	SegTree() noexcept = default;
	SegTree(const SegTree &) = default;
	SegTree(SegTree &&) noexcept = default;
	SegTree &operator=(const SegTree &) = default;
	SegTree &operator=(SegTree &&) noexcept = default;
	
	SegTree(size_t length) :
		m_length(length),
		m_nleaves(calc_nleaves(length)),
		m_data(2 * m_nleaves)
	{}
	
	SegTree(size_t length, const Val &value) : SegTree(length) {
		for (size_t i = 0; i < m_length; ++i) {
			m_data[m_nleaves + i].val = value;
		}
		for (size_t i = m_nleaves - 1; i > 0; --i) {
			m_data[i].val = join(m_data[2*i+0].val, m_data[2*i+1].val);
		}
	}
	
	// std::vector has this constructor only participate in overload resolution if Iter satisfies LegacyInputIterator.
	// To stay compatible with C++17 we can't use concepts, and I'm not gonna write my own sfinae just to check named reqs.
	// Thus, you provide from_iter_tag{} as the first argument to use the constructor taking an iterator range.
	template<typename Iter>
	SegTree(from_iter_tag, Iter first, Iter last) : SegTree(std::distance(first, last)) {
		for (size_t i = 0; i < m_length; ++i) {
			m_data[m_nleaves + i].val = *first;
			++first;
		}
		for (size_t i = m_nleaves - 1; i > 0; --i) {
			m_data[i].val = join(m_data[2*i+0].val, m_data[2*i+1].val);
		}
	}
	
	Val query(size_t first, size_t last) const {
		if (last > m_length) last = m_length;
		if (last <= first) return neutral();
		size_t p = m_nleaves + first;
		size_t q = m_nleaves + last - 1;
		
		Val left = m_data[p].modval();
		
		if (p != q) {
			Val right = m_data[q].modval();
			while (p / 2 != q / 2) {
				if (p % 2 == 0) left = join(left, m_data[p+1].modval());
				if (q % 2 == 1) right = join(m_data[q-1].modval(), right);
				p /= 2; q /= 2;
				left = apply(m_data[p].mod, left);
				right = apply(m_data[q].mod, right);
			}
			left = join(left, right);
		}
		
		while (p /= 2) {
			left = apply(m_data[p].mod, left);
		}
		
		return left;
	}
	
	void update(size_t first, size_t last, const Mod &f) {
		if (last > m_length) last = m_length;
		if (last <= first) return;
		size_t p = m_nleaves + first;
		size_t q = m_nleaves + last - 1;
		
		if (p == q) {
			// k = min k such that (lca >> k) == 0. (aka # of nodes on the path from lca to root)
			size_t k = (std::numeric_limits<unsigned long long>::digits - __builtin_clzll(p));
			
			// lca, p, and q are all the same (leaf) node. Don't down-propagate from leaves.
			for (size_t i = k - 1; i > 0; --i) {
				downprop(p >> i, false);
			}
		} else {
			// m = min m such that (p >> m) == (q >> m). (aka (p^q) >> m == 0) (aka # of nodes on the path from p/2 to lca(p,q))
			size_t m = (std::numeric_limits<unsigned long long>::digits - __builtin_clzll(p ^ q));
			
			size_t lca = p >> m;
			
			// k = min k such that (lca >> k) == 0. (aka # of nodes on the path from lca to root)
			size_t k = (std::numeric_limits<unsigned long long>::digits - __builtin_clzll(lca));
			
			for (size_t i = k - 1; i < k; --i) {
				downprop(lca >> i, false);
			}
			
			for (size_t i = m - 1; i > 0; --i) {
				downprop(p >> i, false);
				downprop(q >> i, false);
			}
		}
		
		
		m_data[p].mod = compose(f, m_data[p].mod);
		
		if (p != q) {
			m_data[q].mod = compose(f, m_data[q].mod);
			while (p / 2 != q / 2) {
				if (p % 2 == 0) m_data[p+1].mod = compose(f, m_data[p+1].mod);
				if (q % 2 == 1) m_data[q-1].mod = compose(f, m_data[q-1].mod);
				p /= 2; q /= 2;
				recalc(p);
				recalc(q);
			}
		}
		
		while (p /= 2) {
			recalc(p);
		}
	}
};
