#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

// Integer 1:
// Implement a signed big integer class that only needs to support simple addition and subtraction

// Integer 2:
// Implement a signed big integer class that supports addition, subtraction, multiplication, and division, and overload related operators

// Do not use any header files other than the following
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

// Do not use "using namespace std;"

namespace sjtu {
class int2048 {
  static const int base = 10000;
  static const int base_digits = 4;
  std::vector<int> digits;
  int sign;

  void normalize();
  static int compare_abs(const int2048 &, const int2048 &);
  static int2048 multiply_abs(const int2048 &, const int2048 &);
  static int2048 divide_abs(const int2048 &, const int2048 &, int2048 *);
  void multiply_small(int);
  void add_small(int);
public:
  // Constructors
  int2048();
  int2048(long long);
  int2048(const std::string &);
  int2048(const int2048 &);

  // The parameter types of the following functions are for reference only, you can choose to use constant references or not
  // If needed, you can add other required functions yourself
  // ===================================
  // Integer1
  // ===================================

  // Read a big integer
  void read(const std::string &);
  // Output the stored big integer, no need for newline
  void print();

  // Add a big integer
  int2048 &add(const int2048 &);
  // Return the sum of two big integers
  friend int2048 add(int2048, const int2048 &);

  // Subtract a big integer
  int2048 &minus(const int2048 &);
  // Return the difference of two big integers
  friend int2048 minus(int2048, const int2048 &);

  // ===================================
  // Integer2
  // ===================================

  int2048 operator+() const;
  int2048 operator-() const;

  int2048 &operator=(const int2048 &);

  int2048 &operator+=(const int2048 &);
  friend int2048 operator+(int2048, const int2048 &);

  int2048 &operator-=(const int2048 &);
  friend int2048 operator-(int2048, const int2048 &);

  int2048 &operator*=(const int2048 &);
  friend int2048 operator*(int2048, const int2048 &);

  int2048 &operator/=(const int2048 &);
  friend int2048 operator/(int2048, const int2048 &);

  int2048 &operator%=(const int2048 &);
  friend int2048 operator%(int2048, const int2048 &);

  friend std::istream &operator>>(std::istream &, int2048 &);
  friend std::ostream &operator<<(std::ostream &, const int2048 &);

  friend bool operator==(const int2048 &, const int2048 &);
  friend bool operator!=(const int2048 &, const int2048 &);
  friend bool operator<(const int2048 &, const int2048 &);
  friend bool operator>(const int2048 &, const int2048 &);
  friend bool operator<=(const int2048 &, const int2048 &);
  friend bool operator>=(const int2048 &, const int2048 &);
};
} // namespace sjtu

#endif

namespace sjtu {

namespace {
const double PI = 3.141592653589793238462643383279502884;

void fft(std::vector<std::complex<double> > &a, bool invert) {
  const int n = static_cast<int>(a.size());
  for (int i = 1, j = 0; i < n; ++i) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) { std::complex<double> t = a[i]; a[i] = a[j]; a[j] = t; }
  }
  for (int len = 2; len <= n; len <<= 1) {
    const double angle = 2.0 * PI / len * (invert ? -1.0 : 1.0);
    const std::complex<double> wlen(std::cos(angle), std::sin(angle));
    for (int i = 0; i < n; i += len) {
      std::complex<double> w(1.0, 0.0);
      for (int j = 0; j < len / 2; ++j) {
        const std::complex<double> u = a[i + j];
        const std::complex<double> v = a[i + j + len / 2] * w;
        a[i + j] = u + v;
        a[i + j + len / 2] = u - v;
        w *= wlen;
      }
    }
  }
  if (invert) for (int i = 0; i < n; ++i) a[i] /= n;
}
}

int2048::int2048() : sign(1) {}
int2048::int2048(long long value) : sign(value < 0 ? -1 : 1) {
  unsigned long long x = value < 0 ? 0ULL - static_cast<unsigned long long>(value) : static_cast<unsigned long long>(value);
  while (x) { digits.push_back(static_cast<int>(x % base)); x /= base; }
}
int2048::int2048(const std::string &s) { read(s); }
int2048::int2048(const int2048 &other) : digits(other.digits), sign(other.sign) {}

void int2048::normalize() {
  while (!digits.empty() && digits.back() == 0) digits.pop_back();
  if (digits.empty()) sign = 1;
}

void int2048::read(const std::string &s) {
  digits.clear(); sign = 1;
  std::size_t begin = 0;
  if (!s.empty() && (s[0] == '-' || s[0] == '+')) { sign = s[0] == '-' ? -1 : 1; begin = 1; }
  while (begin < s.size() && s[begin] == '0') ++begin;
  for (std::size_t end = s.size(); end > begin;) {
    std::size_t start = end >= begin + base_digits ? end - base_digits : begin;
    int x = 0;
    for (std::size_t i = start; i < end; ++i) x = x * 10 + s[i] - '0';
    digits.push_back(x);
    end = start;
  }
  normalize();
}

void int2048::print() { std::cout << *this; }

int int2048::compare_abs(const int2048 &a, const int2048 &b) {
  if (a.digits.size() != b.digits.size()) return a.digits.size() < b.digits.size() ? -1 : 1;
  for (std::size_t i = a.digits.size(); i-- > 0;) if (a.digits[i] != b.digits[i]) return a.digits[i] < b.digits[i] ? -1 : 1;
  return 0;
}

int2048 &int2048::add(const int2048 &rhs) { return *this += rhs; }
int2048 add(int2048 lhs, const int2048 &rhs) { return lhs += rhs; }
int2048 &int2048::minus(const int2048 &rhs) { return *this -= rhs; }
int2048 minus(int2048 lhs, const int2048 &rhs) { return lhs -= rhs; }

int2048 int2048::operator+() const { return *this; }
int2048 int2048::operator-() const { int2048 result(*this); if (!result.digits.empty()) result.sign = -result.sign; return result; }
int2048 &int2048::operator=(const int2048 &rhs) { digits = rhs.digits; sign = rhs.sign; return *this; }

int2048 &int2048::operator+=(const int2048 &rhs) {
  if (rhs.digits.empty()) return *this;
  if (digits.empty()) { *this = rhs; return *this; }
  if (sign == rhs.sign) {
    int carry = 0;
    const std::size_t n = digits.size() > rhs.digits.size() ? digits.size() : rhs.digits.size();
    digits.resize(n);
    for (std::size_t i = 0; i < n || carry; ++i) {
      if (i == digits.size()) digits.push_back(0);
      long long x = digits[i] + carry + (i < rhs.digits.size() ? rhs.digits[i] : 0);
      digits[i] = static_cast<int>(x % base); carry = static_cast<int>(x / base);
    }
  } else {
    int relation = compare_abs(*this, rhs);
    if (relation == 0) { digits.clear(); sign = 1; return *this; }
    const int2048 *large = relation > 0 ? this : &rhs;
    const int2048 *small = relation > 0 ? &rhs : this;
    std::vector<int> result = large->digits;
    int borrow = 0;
    for (std::size_t i = 0; i < result.size(); ++i) {
      int x = result[i] - borrow - (i < small->digits.size() ? small->digits[i] : 0);
      if (x < 0) { x += base; borrow = 1; } else borrow = 0;
      result[i] = x;
    }
    digits.swap(result); sign = large->sign; normalize();
  }
  return *this;
}

int2048 &int2048::operator-=(const int2048 &rhs) { return *this += -rhs; }
int2048 operator+(int2048 lhs, const int2048 &rhs) { return lhs += rhs; }
int2048 operator-(int2048 lhs, const int2048 &rhs) { return lhs -= rhs; }

void int2048::multiply_small(int value) {
  if (value == 0 || digits.empty()) { digits.clear(); sign = 1; return; }
  long long carry = 0;
  for (std::size_t i = 0; i < digits.size() || carry; ++i) {
    if (i == digits.size()) digits.push_back(0);
    long long x = carry + static_cast<long long>(digits[i]) * value;
    digits[i] = static_cast<int>(x % base); carry = x / base;
  }
}
void int2048::add_small(int value) {
  for (std::size_t i = 0; value; ++i) {
    if (i == digits.size()) digits.push_back(0);
    int x = digits[i] + value;
    digits[i] = x % base; value = x / base;
  }
}

int2048 int2048::multiply_abs(const int2048 &a, const int2048 &b) {
  int2048 result;
  if (a.digits.empty() || b.digits.empty()) return result;
  const std::size_t min_size = a.digits.size() < b.digits.size() ? a.digits.size() : b.digits.size();
  if (min_size < 48) {
    result.digits.assign(a.digits.size() + b.digits.size(), 0);
    for (std::size_t i = 0; i < a.digits.size(); ++i) {
      long long carry = 0;
      for (std::size_t j = 0; j < b.digits.size() || carry; ++j) {
        long long x = result.digits[i + j] + carry + (j < b.digits.size() ? static_cast<long long>(a.digits[i]) * b.digits[j] : 0);
        result.digits[i + j] = static_cast<int>(x % base); carry = x / base;
      }
    }
  } else {
    // Split base-10000 limbs into base-100 coefficients before the FFT.  The
    // smaller convolution coefficients keep floating-point rounding safely
    // below one unit even for the largest required products.
    const int convolution_base = 100;
    const std::size_t coefficient_count = 2 * (a.digits.size() + b.digits.size());
    int n = 1;
    while (n < static_cast<int>(coefficient_count)) n <<= 1;
    std::vector<std::complex<double> > fa(n), fb(n);
    for (std::size_t i = 0; i < a.digits.size(); ++i) { fa[2 * i] = a.digits[i] % convolution_base; fa[2 * i + 1] = a.digits[i] / convolution_base; }
    for (std::size_t i = 0; i < b.digits.size(); ++i) { fb[2 * i] = b.digits[i] % convolution_base; fb[2 * i + 1] = b.digits[i] / convolution_base; }
    fft(fa, false); fft(fb, false);
    for (int i = 0; i < n; ++i) fa[i] *= fb[i];
    fft(fa, true);
    std::vector<int> coefficients(coefficient_count + 1);
    long long carry = 0;
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
      long long x = carry + (i < static_cast<std::size_t>(n) ? static_cast<long long>(fa[i].real() + 0.5) : 0);
      coefficients[i] = static_cast<int>(x % convolution_base); carry = x / convolution_base;
    }
    result.digits.resize((coefficients.size() + 1) / 2);
    for (std::size_t i = 0; i < result.digits.size(); ++i) result.digits[i] = coefficients[2 * i] + (2 * i + 1 < coefficients.size() ? convolution_base * coefficients[2 * i + 1] : 0);
  }
  result.normalize(); return result;
}

int2048 &int2048::operator*=(const int2048 &rhs) {
  int result_sign = sign * rhs.sign;
  int2048 result = multiply_abs(*this, rhs);
  if (!result.digits.empty()) result.sign = result_sign;
  *this = result; return *this;
}
int2048 operator*(int2048 lhs, const int2048 &rhs) { return lhs *= rhs; }

int2048 int2048::divide_abs(const int2048 &a, const int2048 &b, int2048 *remainder) {
  if (compare_abs(a, b) < 0) { if (remainder) *remainder = a; return int2048(); }
  if (b.digits.size() == 1) {
    int2048 q; q.digits.resize(a.digits.size()); long long rem = 0;
    for (std::size_t i = a.digits.size(); i-- > 0;) { long long cur = a.digits[i] + rem * base; q.digits[i] = static_cast<int>(cur / b.digits[0]); rem = cur % b.digits[0]; }
    q.normalize(); if (remainder) *remainder = int2048(rem); return q;
  }
  int norm = base / (b.digits.back() + 1);
  int2048 u(a), v(b); u.multiply_small(norm); v.multiply_small(norm);
  u.digits.push_back(0);
  const std::size_t n = v.digits.size(), m = u.digits.size() - n - 1;
  int2048 q; q.digits.assign(m + 1, 0);
  for (std::size_t ii = m + 1; ii-- > 0;) {
    long long numerator = static_cast<long long>(u.digits[ii + n]) * base + u.digits[ii + n - 1];
    long long guess = numerator / v.digits[n - 1];
    long long second = numerator % v.digits[n - 1];
    if (guess >= base) { guess = base - 1; second += v.digits[n - 1]; }
    while (n > 1 && guess * v.digits[n - 2] > static_cast<long long>(base) * second + u.digits[ii + n - 2]) { --guess; second += v.digits[n - 1]; if (second >= base) break; }
    long long borrow = 0, carry = 0;
    for (std::size_t j = 0; j < n; ++j) {
      long long product = guess * v.digits[j] + carry; carry = product / base; product %= base;
      long long x = u.digits[ii + j] - product - borrow;
      if (x < 0) { x += base; borrow = 1; } else borrow = 0;
      u.digits[ii + j] = static_cast<int>(x);
    }
    long long top = u.digits[ii + n] - carry - borrow;
    if (top < 0) {
      --guess; long long c = 0;
      for (std::size_t j = 0; j < n; ++j) { long long x = u.digits[ii + j] + v.digits[j] + c; if (x >= base) { x -= base; c = 1; } else c = 0; u.digits[ii + j] = static_cast<int>(x); }
      u.digits[ii + n] = static_cast<int>(top + base + c);
    } else u.digits[ii + n] = static_cast<int>(top);
    q.digits[ii] = static_cast<int>(guess);
  }
  q.normalize();
  if (remainder) { remainder->digits.assign(u.digits.begin(), u.digits.begin() + n); remainder->sign = 1; long long rem = 0; for (std::size_t i = remainder->digits.size(); i-- > 0;) { long long cur = remainder->digits[i] + rem * base; remainder->digits[i] = static_cast<int>(cur / norm); rem = cur % norm; } remainder->normalize(); }
  return q;
}

int2048 &int2048::operator/=(const int2048 &rhs) {
  int2048 a(*this), b(rhs); a.sign = b.sign = 1;
  int2048 rem; int2048 q = divide_abs(a, b, &rem);
  if (sign != rhs.sign && !rem.digits.empty()) { q.add_small(1); q.sign = -1; } else if (!q.digits.empty()) q.sign = sign * rhs.sign;
  *this = q; return *this;
}
int2048 operator/(int2048 lhs, const int2048 &rhs) { return lhs /= rhs; }

int2048 &int2048::operator%=(const int2048 &rhs) {
  int2048 quotient = *this / rhs;
  *this -= quotient * rhs;
  return *this;
}
int2048 operator%(int2048 lhs, const int2048 &rhs) { return lhs %= rhs; }

std::istream &operator>>(std::istream &in, int2048 &value) { std::string s; in >> s; value.read(s); return in; }
std::ostream &operator<<(std::ostream &out, const int2048 &value) {
  if (value.sign < 0) out << '-';
  if (value.digits.empty()) return out << '0';
  out << value.digits.back();
  for (std::size_t i = value.digits.size() - 1; i-- > 0;) { char buf[5]; std::snprintf(buf, sizeof(buf), "%04d", value.digits[i]); out << buf; }
  return out;
}

bool operator==(const int2048 &a, const int2048 &b) { return a.sign == b.sign && a.digits == b.digits; }
bool operator!=(const int2048 &a, const int2048 &b) { return !(a == b); }
bool operator<(const int2048 &a, const int2048 &b) { if (a.sign != b.sign) return a.sign < b.sign; int x = int2048::compare_abs(a, b); return a.sign > 0 ? x < 0 : x > 0; }
bool operator>(const int2048 &a, const int2048 &b) { return b < a; }
bool operator<=(const int2048 &a, const int2048 &b) { return !(b < a); }
bool operator>=(const int2048 &a, const int2048 &b) { return !(a < b); }

} // namespace sjtu
