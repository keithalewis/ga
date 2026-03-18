// ga.h - Grassmann algebra
#pragma once

#include <bitset>
#include <charconv>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace fms::ga
{
	// Product of points P_i, 0 <= i < N, represented as a bitmask of non-zero P_i.
	template <std::size_t N = 64>
	using blade = std::bitset<N>;

	// Index of largest set bit
	template<std::size_t N = 64>
	constexpr std::size_t order(const blade<N>& p)
	{
		if (p.none()) {
			return 0; // no bits set
		}

		for (std::size_t i = N; i-- > 0;) {
			if (p.test(i)) {
				return i + 1;
			}
		}

		return 0; // no bits set
	}
	static_assert(order(blade<5>(0b0)) == 0);
	static_assert(order(blade<5>(0b1)) == 1);
	static_assert(order(blade<5>(0b10)) == 2);
	static_assert(order(blade<5>(0b10000)) == 5);
	static_assert(order(blade<5>(0b10100)) == 5);
	static_assert(order(blade<5>(0b10101)) == 5);

	// P_i
	template<std::size_t N = 64>
	constexpr blade<N> P(std::size_t i)
	{
		return blade<N>().set(i);
	}
	static_assert(order(P(0)) == 1);
	static_assert(order(P(1)) == 2);
	static_assert(order(P(4)) == 5);

	// 0b1101 -> "0.2.3"
	// non-zero blade bit positions
	template <std::size_t N = 64>
	constexpr std::to_chars_result to_chars(std::string_view s, blade<N> p, char sep = '.')
	{
		char* b = const_cast<char*>(s.data());
		char* e = b + s.size();
		char dot = 0;
		for (std::size_t i = 0; i < N; ++i) {
			if (p.test(i)) {
				if (dot) {
					*b++ = dot;
				}
				auto [ptr, ec] = std::to_chars(b, e, i);
				if (ec != std::errc()) {
					return { b, ec }; // error
				}
				b = ptr;
				dot = sep;
			}
		}

		return { b, std::errc() };
	}
	//constexpr auto xxx = to_chars("0.2.3");
	/*
	namespace {
		constexpr char buf[17];
		constexpr std::string_view v(buf, buf + 17);
		constexpr auto s = to_chars(std::string_view(buf), blade<5>(0b1101));
		static_assert(s.ec == std::errc());
		static_assert(s.ptr == s.data() + 5);
		static_assert(std::string_view(s.data(), s.ptr - s.data()) == "0.2.3");
	}
	*/

	// "0.2.3" -> 0b1101 
	// non-zero blade bit positions
	template <std::size_t N = 64>
	constexpr blade<N> from_string(std::string_view s, char sep = '.')
	{
		blade<N> p;
		std::size_t start = 0;
		while (start < s.size()) {
			auto end = s.find(sep, start);
			if (end == std::string_view::npos) {
				end = s.size();
			}
			std::size_t index;
			auto [ptr, ec] = std::from_chars(s.data() + start, s.data() + end, index);
			if (ec == std::errc() && index < N) {
				p.set(index);
			}
			else {
				return blade<N>(); // invalid input
			}
			start = end + 1;
		}

		return p;
	}
	static_assert(from_string<5>("") == blade<5>());
	static_assert(from_string<5>("0") == blade<5>(0b1));
	static_assert(from_string<5>("1") == blade<5>(0b10));
	static_assert(from_string<5>("0.2.3") == blade<5>(0b1101));

	// Number of non-zero bits.
	template <std::size_t N = 64>
	constexpr std::size_t grade(const blade<N>& p)
	{
		return p.count();
	}
	static_assert(grade(blade<5>()) == 0);
	static_assert(grade(blade<5>("10101")) == 3);
	static_assert(grade(blade<5>(0b11111)) == 5);

	// Grassmann exterior product including sign.
	template <std::size_t N>
	constexpr std::pair<int, blade<N>> exterior_product(const blade<N>& p, const blade<N>& q)
	{
		if ((p & q).any()) // PP = 0
		{
			return { 0, blade<N>(0) }; // zero element
		}

		int sign = 0;
		int seen = 0;

		for (std::size_t i = 0; i < N; ++i) {
			if (p[i]) {
				++seen;
			}
			if (q[i]) {
				sign += seen;
			}
		}

		return { sign % 2 ? 1 : -1, p | q };
	}
	// P0 P0 = 0
	static_assert(exterior_product(blade<5>("1"), blade<5>("1"))
		== std::make_pair(0, blade<5>("")));
	// P0 P1 = P0 P1
	static_assert(exterior_product(blade<5>("1"), blade<5>("10"))
		== std::make_pair(1, blade<5>("11")));
	// P1 P0 = - P0 P1
	static_assert(exterior_product(blade<5>("10"), blade<5>("1"))
		== std::make_pair(-1, blade<5>("11")));
	// P1 P3 v P0 P2 = P0 P1 P2 P3
	static_assert(exterior_product(blade<5>("1010"), blade<5>("0101"))
		== std::make_pair(1, blade<5>("1111")));

	// Comparison operator for blades.
	template <std::size_t N>
	struct blade_less {
		static constexpr bool operator()(const blade<N>& a, const blade<N>& b) {
			//return a.to_string() < b.to_string();
			return a.to_ullong() < b.to_ullong(); // only if N <= 64
		}
	};
	static_assert(blade_less<5>()(blade<5>("1"), blade<5>("10")) == true);
	static_assert(blade_less<5>()(blade<5>("10"), blade<5>("1")) == false);

	// General element of grassmann algebra.
	template <std::size_t N = 64, typename T = double>
	//requires std::is_floating_point_v<T>
	class extent {
		std::map<blade<N>, T, blade_less<N>> x_;
		// remove blades with zero coefficients
		constexpr extent& trim()
		{
			for (auto it = x_.begin(); it != x_.end();) {
				if (it->second == T(0)) {
					it = x_.erase(it);
				}
				else {
					++it;
				}
			}

			return *this;
		}
	public:
		constexpr extent() = default;
		constexpr extent(const extent&) = default;
		constexpr extent& operator=(const extent&) = default;
		constexpr extent(extent&&) = default;
		constexpr extent& operator=(extent&&) = default;
		constexpr ~extent() = default;

		// scalar
		constexpr extent(const T& x)
		{
			x_[blade<N>()] = x;
		}

		// x P_i
		constexpr extent(const T& x, std::size_t i)
		{
			x_[P(i)] = x;
		}

		// x P
		constexpr extent(const T& x, const blade<N>& p)
		{
			x_[p] = x;
		}

		// Point corresponding to sum_i x_i P_i with weight sum_i x_i.
		template<std::size_t N_ = N>
			requires (N_ <= N || N_ == std::dynamic_extent)
		constexpr extent(const std::span<const T, N_>& x)
		{
			if (x.size() > N) {
				throw std::invalid_argument("x must have size less than or equal to N");
			}

			for (std::size_t i = 0; i < x.size(); ++i) {
				if (x[i]) {
					operator+=(extent(x[i], i));
				}
			}
		}

		constexpr std::size_t size() const
		{
			return x_.size();
		}
		constexpr bool empty() const
		{
			return x_.empty();
		}
		constexpr auto begin() const
		{
			return x_.begin();
		}
		constexpr auto end() const
		{
			return x_.end();
		}

		// access coefficient by basis blade
		constexpr T& operator[](const blade<N>& p) {
			return x_[p];
		}
		constexpr const T& operator[](const blade<N>& p) const {
			return x_.at(p);
		}

		constexpr bool operator==(const extent& b) const
		{
			return x_ == b.x_;
		}

		// {(x, 0)}
		constexpr bool is_scalar() const
		{
			return size() == 1 and x_.begin()->first.none();
		}
		// {(x1, b1), (x2, b10), (x3, b100) ...}
		constexpr bool is_point() const
		{
			bool point = true;

			for (const auto& [p, x] : x_) {
				if (p.count() != 1) {
					point = false;
					break;
				}
			}

			return point;
		}
		// Maximum dimension of non-zero blades.
		constexpr std::size_t depth() const
		{
			std::size_t d = 0;

			for (const auto& [p, x] : x_) {
				if (x) {
					d = (std::max)(d, order(p));
				}
			}

			return d;
		}
		// total weight of each grade
		constexpr std::array<T, N + 1> weight() const
		{
			std::array<T, N + 1> w{};

			for (const auto& [p, x] : x_) {
				w[p.count()] += x;
			}

			return w;
		}

		// complement
		constexpr extent& operator~()
		{
			extent _x;

			for (const auto& [p, x] : x_) {
				_x[~p] = x;
			}
			_x.x_.swap(x_);

			return *this;
		}
		// additive inverse
		constexpr extent& operator-()
		{
			for (auto& [p, x] : x_) {
				x_[p] = -x;
			}

			return *this;
		}
		constexpr extent& operator+=(const extent& b)
		{
			for (const auto& [p, x] : b) {
				x_[p] += x;
			}

			return trim();
		}
		// Add x(P_i - P_0)
		constexpr extent& add(T x, std::size_t i)
		{
			operator+=(extent(x, i) - extent(x, 0));

			return *this;
		}
		constexpr extent& operator-=(const extent& b)
		{
			for (const auto& [p, x] : b) {
				x_[p] -= x;
			}

			return trim();
		}
		// exterior product
		constexpr extent& operator|(const extent& e)
		{
			extent _x;

			// distributive law
			for (const auto& [p, x] : x_) {
				for (const auto& [q, y] : e) {
					const auto [sign, pq] = exterior_product(p, q);
					T xy = x * y;
					if (pq.any() and xy) { // non-zero
						_x[pq] += sign * xy;
					}
				}
			}
			std::swap(_x, x_);

			return trim();
		}
		constexpr extent& operator*=(const T& c)
		{
			for (auto& [p, x] : x_) {
				x *= c;
			}

			return trim();
		}

		// Has single blade same as e.
		constexpr blade<N> congruent(const extent& e) const
		{
			const auto b = begin()->first;
			const auto eb = e.begin()->first;

			return (size() == 1 and e.size() == 1 and b == eb)
				? b : blade<N>{};
		}
		// quotient of coefficients if congruent, otherwise NaN
		constexpr T operator/(const extent& e) const
		{
			const blade<N> f = congruent(e);

			return x_[f] / e.x_[f];
		}
		static_assert(extent(1, blade()) / extent(2, blade()) == 1 / 2);
		static_assert(extent(1, blade(0b1)) / extent(2, blade(0b1)) == 1 / 2);

		// Return coefficients in terms of P_i
		constexpr extent<N, T> fit(std::span<const extent<N, T>> P) const
		{
			extent<N, T> Q0; // P0 ... P_{N-1}
			std::vector<extent<N, T>> Q(P.size());

			for (std::size_t i = 0; i < P.size(); ++i) {
				for (std::size_t j = 0; j < P.size(); ++j) {
					Q0 |= P[j];
					Q[i] |= (i == j) ? *this : P[j];
				}
			}

			fms::ga::extent<N, T> R;
			for (std::size_t i = 0; i < P.size(); ++i) {
				R += fms::ga::extent(Q[i] / Q0, i);
			}
			// assert *this == sum_i R[i] P[i]

			return R;
		}

	};
	// empty vs scalar 0.
	//static_assert(extent<>() != extent<>(0));

	template <std::size_t N = 64, typename T = double>
	constexpr extent<N, T> _0 = extent<N, T>(T(0));
}

// complement
template <std::size_t N, typename T = double>
	requires std::is_floating_point_v<T>
constexpr fms::ga::extent<N, T> operator~(fms::ga::extent<N, T> a)
{
	return ~a;
}

// exterior product
template <std::size_t N, typename T = double>
	requires std::is_floating_point_v<T>
constexpr fms::ga::extent<N, T> operator|(fms::ga::extent<N, T> a, const fms::ga::extent<N, T>& b)
{
	return a |= b;
}

// interior product
template <std::size_t N, typename T = double>
	requires std::is_floating_point_v<T>
constexpr fms::ga::extent<N, T> operator&(fms::ga::extent<N, T> a, const fms::ga::extent<N, T>& b)
{
	return ~((~a) | (~b));
}

template <std::size_t N, typename T = double>
	requires std::is_floating_point_v<T>
constexpr fms::ga::extent<N, T> operator+(fms::ga::extent<N, T> a, const fms::ga::extent<N, T>& b)
{
	return a += b;
}

template <std::size_t N, typename T = double>
	requires std::is_floating_point_v<T>
constexpr fms::ga::extent<N, T> operator-(fms::ga::extent<N, T> a, const fms::ga::extent<N, T>& b)
{
	return a -= b;
}

// scalar multiplication
template <std::size_t N, typename T = double>
	requires std::is_floating_point_v<T>
constexpr fms::ga::extent<N, T> operator*(fms::ga::extent<N, T> a, T b)
{
	return a *= b;
}

template <std::size_t N, typename T = double>
	requires std::is_floating_point_v<T>
constexpr fms::ga::extent<N, T> operator*(T b, fms::ga::extent<N, T> a)
{
	return a * b;
}

// scalar division
template <std::size_t N, typename T = double>
	requires std::is_floating_point_v<T>
constexpr fms::ga::extent<N, T> operator/(fms::ga::extent<N, T> a, T b)
{
	return a /= b;
}

// congruent division
template <std::size_t N, typename T = double>
	requires std::is_floating_point_v<T>
constexpr fms::ga::extent<N, T> operator/(fms::ga::extent<N, T> a, fms::ga::extent<N, T> b)
{
	return a /= b;
}

