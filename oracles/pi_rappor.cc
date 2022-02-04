#include "pi_rappor.h"
#include "util/util.h"

using namespace std;

vector<int> PIAsymmetricRAPPOR::num_to_vec(int x, int l) {
  vector<int> ret(l);
  for (int i = 0; i < l; ++i) {
    ret[l-1-i] = x % q;
    x /= q;
  }
  return ret;
}

int PIAsymmetricRAPPOR::vec_to_num(const vector<int> &v) {
  int x = 0, qpow = 1;
  for (int i = 0; i < v.size(); ++i) {
    x += v[v.size()-1-i]*qpow;
    qpow *= q;
  }
  return x;
}

PIAsymmetricRAPPOR::PIAsymmetricRAPPOR(int K_, double eps, bool debug, uint32_t seed)
  : PrivateFrequencyOracle(K_, eps, debug, seed) {
  p1 = 1.0 / (exp(eps) + 1);
  p2 = 0.5;
  // to get DP, need 1/q >= (1 - p2) / (1 + p2(exp(eps) - 1))
  // i.e. q <= (1 + p2(exp(eps) - 1)) / (1 - p2)
  q = (1 + p2*(exp(eps) - 1)) / (1.0 - p2);
  if (q%2==0 && q!=2)
    --q;
  while (!Util::is_prime(q))
    q -= 2;
  zero_threshold = q-1; // value >= zero_threshold tells us to flip 0->1; else don't flip
  p1 = 1 - 1.0 * zero_threshold / q;
  int qpow = 1;
  t = 0;
  qpows.push_back(qpow);
  while (qpow <= K) {
    ++t;
    qpow *= q;
    qpows.push_back(qpow);
  }
  K = qpow - 1;
  qinv = Util::compute_inverse_table(q, rng);
  alpha = 1.0 / (1.0 - p1 - p2);
  beta = -p1 / (1.0 - p1 - p2);
  unif_intq = boost::uniform_int<>(0, q-1);
  unif_int0 = boost::uniform_int<>(0, zero_threshold - 1);
  unif_int1 = boost::uniform_int<>(zero_threshold, q - 1);
  if (debug) cerr << "PIAsymetricRAPPOR zero_threshold is " << zero_threshold << endl;
  if (debug) cerr << "PIAsymetricRAPPOR alpha,beta=" << alpha << "," << beta << endl;
  if (debug) cerr << "PIAsymetricRAPPOR exp(eps)+1=" << exp(eps) + 1 << endl;
  if (debug) cerr << "PIAsymetricRAPPOR K,q,t=" << K << "," << q << "," << t << endl;
}

Message PIAsymmetricRAPPOR::local_randomizer(int x) {
  vector<int> ret(t + 1, 0), v = num_to_vec(x + 1, t);
  int dot_prod = 0;
  for (int i = 0; i < t; ++i) {
    ret[i] = unif_intq(rng);
    (dot_prod += (int64_t)ret[i]*v[i]) %= q;
  }
  if (unif_fraction(rng) <= p2) 
    ret[t] = (q + unif_int0(rng) - dot_prod) % q;
  else 
    ret[t] = (q + unif_int1(rng) - dot_prod) % q;
  return Message(vec_to_num(ret));
}

double PIAsymmetricRAPPOR::estimate_freq(int x, const vector<Message> &messages) {
  vector<int> v = num_to_vec(x + 1, t);
  int cnt = 0;
  for (Message m : messages) {
    int dot_prod = 0;
    vector<int> u = num_to_vec(m.read(), t+1);
    for (int i = 0; i < v.size(); ++i)
      (dot_prod += (int64_t)u[i]*v[i]) %= q;
    (dot_prod += u[t]) %= q;
    if (dot_prod >= zero_threshold)
      ++cnt;
  }
  return alpha*cnt + beta*messages.size();
}

// y[u] is the number of messages received equal to u
vector<int> PIAsymmetricRAPPOR::dp_bottom_up(const vector<int> &y) {
  // dynamic programming
  // f(a,b,z) is sum_{u=(u',r) : pref(u')=a, <suff(u'),b> + r = z} y_u
  int N = qpows[t]*q; // mem table requirement for any fixed l = length(a)
  vector<int> last(N), next(N);
  
  // answers for base cases when l==t
  for (int a = 0; a < qpows[t]; ++a)
    for (int z = 0; z < q; ++z)
      last[a*q + z] = y[a*q + z];
  
  int prevA = qpows[t], prevB = 1, curA, curB;
  for (int length = t-1; length >= 0; --length) {
    fill(next.begin(), next.end(), 0);
    curA = qpows[length], curB = qpows[t - length];
    for (int b = 0; b < curB; ++b) {
      int val = b % qpows[t-length-1];
      int first_b_digit = (b - val) / qpows[t-length-1];
      assert(first_b_digit*qpows[t-length-1]+val == b);
      for (int a = 0; a < curA; ++a)
	for (int z = 0; z < q; ++z)
	  for (int64_t d = 0; d < q; ++d) 
	    next[b*curA*q + a*q + z] += last[val*prevA*q + (a*q+d)*q + (q*d + z - first_b_digit*d) % q];
    }
    swap(last, next);
    prevA = curA;
    prevB = curB;
  }
  vector<int> ret(K);
  for (int i = 1; i <= K; ++i)
    for (int z = zero_threshold; z < q; ++z)
      ret[i - 1] += last[i*q + z];
  return ret;
}

vector<double> PIAsymmetricRAPPOR::estimate_all_freqs(const vector<Message> &messages) {
  vector<double> ret(K);
  if (messages.size() < (int64_t)q*q*q)  // naive method faster if Knt < Kq^3 t
    for (int x = 0; x < K; ++x)
      ret[x] = estimate_freq(x, messages);
  else {
    vector<int> y(qpows[t]*q);
    for (Message m : messages) 
      y[m.read()]++;
    vector<int> T = dp_bottom_up(y);
    for (int x = 0; x < K; ++x) 
      ret[x] = alpha*T[x] + beta*messages.size();
  }
  return ret;
}
