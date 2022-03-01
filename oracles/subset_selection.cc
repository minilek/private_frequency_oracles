#include "subset_selection.h"
#include <bits/stdc++.h>

using namespace std;

SubsetSelection::SubsetSelection(int K_, double eps, bool debug_, uint32_t seed) 
  : PrivateFrequencyOracle(K_, eps, debug_, seed) {
  //    K = max(d, K);
  d = (int)ceil(1.0*K / (exp(eps) + 1));
  unif_int = boost::uniform_int<>(0, K-1);
  
  // probability of sending a message in Z_{k,d}^i
  p = exp(eps)*d / (exp(eps)*d + K - d);
  
  // for y!=x, what's the probability that Message(y) has a 1 in the xth bit?
  // it is q = p*(d-1)/(K-1) + (1-p)*d/(K-1)
  double q = p*(d-1)/(K-1) + (1-p)*d/(K-1);
  
  // p = Prob(send vector in Z_{k,d}^i) = (exp(eps)*d / (exp(eps)*d + K - d))
  //
  // unbiased estimator constraints for alpha*[[y in Z_{k,d}^i]] * beta:
  // alpha*p + beta = 1 
  // alpha*q + beta = 0
  //
  // ==> beta = -alpha*q
  // ==> beta = 1 - p*alpha
  // ==> 0 = 1 + (q-p)*alpha, so alpha = 1/(p-q)
  alpha = 1. / (p - q);
  beta = -alpha * q;
  if (debug) cerr << "SS alpha,beta,d=" << alpha << "," << beta << "," << d << endl;
}

// encoding of message in O(d log k) bits. O(d log(k/d)) is possible with O(d) time
// encoding and decoding of messages, but we don't make that optimization here since we
// are not comparing encoding lengths empirically; for the optimized method, see
//
// Larry Carter, Robert W. Floyd, John Gill, George Markowsky, Mark N. Wegman:
// Exact and Approximate Membership Testers. STOC 1978: 59-65,
// the section entitled "Exact Membership Tester 2" (thanks to Huacheng Yu for pointing it out)
// We sketch below how it works; will assume d divides k (not needed, but eases exposition)
// Summary: for subset S of size d, equipartition the universe into d groups of size k/d
// then assign each x in S to the group it's in. Sort the x's according to their
// group number floor(x/(k/d)), which can be done in O(d) time using Counting Sort. Create a
// bitvector b[1..2d-1] of length d intialized to all 0s, and set j, t = 1. Also
// create an array A[1..d], whose elts will always be in {0,...,k/d-1}. Then do:
//
// for i=1..2d-1:
//   let x be the j'th x in S, in sorted order by group number
//   if x is in the t'th group: b[i] = 1, A[j] = x%(k/d), ++j
//   else: b[i] = 0, ++t
// return (A,b) as the encoding of S
// 
// encoding length is d*log(k/d) + 2d - 1 bits.
Message SubsetSelection::local_randomizer(int x) {
  vector<int> ret;
  unordered_set<int> s;
  if (unif_fraction(rng) <= p)
    s.insert(x);
  while (s.size() < d) {
    int y;
    while (true) {
      y = unif_int(rng);
      if (s.find(y) == s.end())
	break;
    }
    s.insert(y);
  }
  for (unordered_set<int>::iterator iter = s.begin(); iter != s.end(); ++iter)
    ret.push_back(*iter);
  return Message(ret);
}

double SubsetSelection::estimate_freq(int x, const vector<Message> &messages) {
  int cnt = 0;
  for (Message m : messages)
    for (int i = 0; i < d; ++i)
      if (m[i] == x) {
	++cnt;
	break;
      }
  return alpha*cnt + beta*messages.size();
}

vector<double> SubsetSelection::estimate_all_freqs(const vector<Message> &messages) {
  vector<double> cnts(K, 0);
  for (Message m : messages)
    for (int i = 0; i < d; ++i)
      cnts[m[i]]++;
  for (int x = 0; x < K; ++x)
    cnts[x] = alpha*cnts[x] + beta*messages.size();
  return cnts;
}
