#include <bits/stdc++.h>
#include <boost/random.hpp>
#include "util.h"

using namespace std;

void Util::generate_divisors(const vector< pair<int,int> > &prime_factorization, vector<int> &ret, int at, int sofar) {
  if (at == prime_factorization.size()) {
    ret.push_back(sofar);
  } else {
    generate_divisors(prime_factorization, ret, at+1, sofar);
    for (int i = 0; i < prime_factorization[at].second; ++i) {
      sofar *= prime_factorization[at].first;
      generate_divisors(prime_factorization, ret, at+1, sofar);
    }
  }
}
  
// calculate generator for multiplicative group of F_q, q prime
int Util::primitive_root(int q, boost::random::mt19937 &rng) {
  int val = q-1;
  vector< pair<int, int> > prime_factorization; // factorize q-1;
  for (int d = 2; d*d <= q; ++d)
    if (val % d == 0) {
      int cnt = 0;
      while (val % d == 0) {
	++cnt;
	val /= d;
      }
      prime_factorization.push_back(make_pair(d, cnt));    
    }
  if (val != 1)
    prime_factorization.push_back(make_pair(val, 1));
  vector<int> divisors;
  generate_divisors(prime_factorization, divisors, 0, 1);
  boost::uniform_int<> unif_int = boost::uniform_int<>(1, q - 1);
  while (true) {
    // a random element of F_p^* has probability 1-o(1) of being a generator
    int z = unif_int(rng);
    bool bad = false;
    for (int div : divisors) {
      if ((div == 1) || (div == q-1))
	continue;
      if (fast_mod_pow(z, div, q) == 1) {
	bad = true;
	break;
      }
    }
    if (!bad) 
      return z;
  }
}  


vector<int> Util::hadamard_transform(vector<int> x) {
  int n = x.size();
  if (n == 1)
    return x;
  else {
    assert((n&1) == 0);
    vector<int> ret(n);
    vector<int> xl(n/2);
    vector<int> xh(n/2);
    for (int i = 0; i < n/2; ++i)
      xl[i] = x[i];
    for (int i = n/2; i < n; ++i)
      xh[i-n/2] = x[i];
    xl = hadamard_transform(xl);
    xh = hadamard_transform(xh);
    for (int i = 0; i < n/2; ++i) {
      ret[i] = xl[i] + xh[i];
      ret[n/2+i] = xl[i] - xh[i];
    }
    return ret;
  }
}

int Util::fast_mod_pow(int x, int n, int q) {
  if (n==0)
    return 1;
  else {
    int z = fast_mod_pow(x, n/2, q);
    z = ((int64_t)z * z) % q;
    if (n&1 == 1)
      z = ((int64_t)z * x) % q;
    return z;
  }
}   

vector<int> Util::compute_inverse_table(int q, boost::random::mt19937 &rng) {
  // primitive root of finite field F_q, allows us to calculate table of all
  // compute table of all F_q multiplicative inverses in O(q) time; naively
  // would be O(q log q), which is O(k log k) for t=2 since q=Theta(k) for t=2
  assert(is_prime(q));
  int z = primitive_root(q, rng);
  vector<int> zpows(q), ret(q);
  zpows[0] = 1;
  for (int i = 1; i < q; ++i) 
    zpows[i] = ((int64_t)zpows[i-1] * z) % q;
  for (int i = 0; i < q-1; ++i) 
    ret[zpows[i]] = zpows[q-i-1];
  return ret;
}

int Util::first_non_zero(const vector<int> &v) {
  int i = 0;
  while ((i<v.size()) && !v[i])
    ++i;
  if (i == v.size())
    return -1;
  else
    return i;
}

// normal a nonzero vector so that it's canonical (first non-entry is a 1)
void Util::canonicalize(vector<int> &v, int q, const vector<int> &qinv) {
  int i = first_non_zero(v);
  assert(i != -1);
  int g = qinv[v[i]];
  for (int j = i; j < v.size(); ++j)
    v[j] = ((int64_t)v[j] * g) % q;
}

vector<int> Util::canonicalized_index_to_vec(int j, int l, int q) {
  // given an index in {0,...,K-1}, return the corresponding canonical vector
  // recall a canonical vector is a vector in {0,...,q-1}^t which is nonzero,
  // and its first nonzero entry (going from smallest vector index to largest)
  // has the value 1
  int tot = 0;
  int qpow = 1;
  int num_zeroes = l - 1;
  while (j >= tot + qpow) {
    tot += qpow;
    qpow *= q;
    num_zeroes -= 1;
  }
  vector<int> ret(l);
  ret[num_zeroes] = 1;
  int r = j - tot;
  int at = l - 1;
  while (r > 0) {
    ret[at] = r % q;
    r /= q;
    at -= 1;
  }
  return ret;
}

int Util::canonicalized_vec_to_index(const vector<int> &v, int q, const vector<int> &qpows) {
  // given a canonical vector in {0,...,q-1}^t, return its index in 
  // {0,...,K-1}
  int leftmost_one = Util::first_non_zero(v);
  assert((leftmost_one!=-1) && (v[leftmost_one]==1));
  int ret = (qpows[v.size() - 1 - leftmost_one] - 1) / (q - 1);
  int c = 0;
  for (int i = leftmost_one + 1; i < v.size(); ++i) 
    c = (c*q) + v[i];
  return ret + c;
}

// vindex is index (1-based indexing) to a canonical vector v in F_q^l
// returns a 3dim-vector (v[0], g, vsuff_index), where g is the value of the
// first non-zero entry of vsuffix (vector v with v[0] removed) or 1 if
// vsuffix is the 0 vector. vsuff_index is the index of vsuffix amongst all
// canonical vectors in F_q^{l-1}
vector<int> Util::decompose_canonical_vector(int vindex, int l, int q,
					     const vector<int> &qpows, const vector<int> &qinv) {
  assert(l > 0);
  vector<int> v(l), vsuffix(l - 1); // vindex and its suffix as vectors
  if (vindex) {
    v = Util::canonicalized_index_to_vec(vindex - 1, l, q);
    assert(Util::canonicalized_vec_to_index(v, q, qpows) == vindex - 1);
    for (int i = 1; i < v.size(); ++i)  
      vsuffix[i - 1] = v[i];
  }
  int g = 1; // value of first nonzero entry in vbsuffix, else 1 if it's all zeroes
  bool vsuffix_all_zeroes = true;
  for (int i = 0; i < vsuffix.size(); ++i)
    if (vsuffix[i]) {
      g = vsuffix[i];
      vsuffix_all_zeroes = false;
      break;
    }
  if (!vsuffix_all_zeroes)
    Util::canonicalize(vsuffix, q, qinv); // multiplies vsuffix through by g^{-1} entry-wise (mod q)
  int vsuff_index = 0; // a number in the range [0, (q^{l-1}-1)/(q-1) + 1]
  if (!vsuffix_all_zeroes) {
    vsuff_index = Util::canonicalized_vec_to_index(vsuffix, q, qpows) + 1;
    assert(Util::canonicalized_index_to_vec(vsuff_index - 1, l - 1, q) == vsuffix);
  }
  assert(vsuff_index <= (qpows[l - 1] - 1) / (q - 1) + 1);
  return vector<int>( { v[0], g, vsuff_index } );
}

int Util::mod_dot_product(const vector<int> &u, const vector<int> &v, int q) {
  assert(u.size() == v.size());
  int ret = 0;
  for (int i = 0; i < u.size(); ++i)
    (ret += ((int64_t)u[i]*v[i]) % q) %= q;
  return ret;
}

bool Util::is_prime(int q) {
  if (q < 2) return false;
  else if (q == 2) return true;
  else if (q % 2 == 0) return false;
  else {
    for (int d = 3; d < q; d += 2) {
      if ((int64_t)d*d > q) break;
      if (q % d == 0) return false;
    }
    return true;
  }
}

// pick a random nonzero vector u in F_q^t amongst those with <u,v> = s
// we assume here that v is not zero (else this isn't possible for all s)
vector<int> Util::randvec_with_specified_dot_product(const vector<int>& v, int s, int q,
						     const vector<int> &qinv, boost::random::mt19937 &rng) {
  boost::uniform_int<> unif_int = boost::uniform_int<>(0, q - 1);
  vector<int> ret(v.size());
  // pick ret to be a random element of (F_q)^t \ {0} s.t. <ret,v> = s mod q
  while (true) {
    for (int i = 0; i < v.size(); ++i)
      ret[i] = unif_int(rng);
    if (accumulate(ret.begin(), ret.end(), 0) == 0)
      continue;
    
    int i = first_non_zero(v);
    assert(i != -1);
    
    // compute the dot product of ret with v, ignoring index i
    int z = 0;
    for (int j = 0; j < v.size(); ++j)
      if (j != i)
	z = (z + (int64_t)v[j]*ret[j]) % q;
    // now set ret[i] to make the dot product s
    // we need z + ret[i]*v[i] = s, so ret[i] = (s - z)*v[i]^{-1}
    ret[i] = ((int64_t)(q + s - z) * qinv[v[i]]) % q;
    // be careful; ret might now be the 0 vector. try again if it is.
    if (accumulate(ret.begin(), ret.end(), 0))
      break;
  }
#ifdef DEBUG
  // make sure the dot product is actually what we said it would be
  int calc = 0;
  for (int i = 0; i < v.size(); ++i)
    (calc += (int64_t)v[i]*ret[i]) %= q;
  assert(calc == s);
#endif
  return ret;
}

