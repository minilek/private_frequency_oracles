#include "hadamard_response.h"
#include "util/util.h"

using namespace std;

HadamardResponse::HadamardResponse(int K_, double eps, bool debug, uint32_t seed) 
  : PrivateFrequencyOracle(K_, eps, debug, seed) {
  int logB = 1;
  double target = min(exp(eps), 2.0*K_);
  while (1 << logB < target)
    logB++;
  logB--;
  B = 1 << logB;
  logb = 0;
  while (1<<logb < K_/B + 1)
    logb++;
  b = 1 << logb;
  K = B * (b-1);
  // p = Prob(send vector not in my set)
  // my set has size b/2
  // total number of elts not in my set is B*b - b/2
  // Prob(pick an element from my set) should be exp(eps) times larger than elt not in my set
  // need p*(b/2)*exp(eps) + p*(B - b/2) = 1
  // solves to p * [(b/2)*exp(eps) + B*b - b/2] = 1
  // so p = 1.0 / ( (b/2)*(exp(eps) - 1 + 2B) )   
  p = 1.0 / ( (b/2.0)*(exp(eps) - 1 + 2*B) );
  
  unif_int = boost::uniform_int<>(0, B-1);
  unif_intb = boost::uniform_int<>(0, b-1);
  unif_intB = boost::uniform_int<>(0, B-2);
  if (debug) cerr << "HadamardResponse K,B,b=" << K << "," << B << "," << b << endl;
}

Message HadamardResponse::local_randomizer(int x) {
  int block = x / (b-1);
  if (unif_fraction(rng) <= p*(b/2)*(exp(eps)+1)) {
    x = (x % (b-1)) + 1;
    int ret = 0, parity = 0;
    for (int i = 0; i < logb; ++i, x >>= 1) 
      if (x == 1)  { // i is the most significant bit
	if (unif_fraction(rng) <= exp(eps) / (exp(eps) + 1))
	  ret += parity * (1<<i), parity = 0;
	else
	  ret += !parity * (1<<i), parity = 1;
      } else if (unif_fraction(rng) <= 0.5)
	parity ^= x&1, ret += 1<<i;
    return Message(block*b + ret);
  } else {
    // otherwise pick a random element not in my block
    int rnd_block = unif_intB(rng);
    if (rnd_block >= block)
      rnd_block++;
    return Message(rnd_block*b + unif_intb(rng));
  }
  
}

double HadamardResponse::estimate_freq(int x, const vector<Message> &messages) {
  int block = x / (b-1), id = (x % (b-1)) + 1;
  int cnt1 = 0, cnt2 = 0;
  for (Message m : messages) {
    int u = m.read();
    if (u/b == block) {
      ++cnt2;
      if (!__builtin_parity(id & (u % b)))
	++cnt1;
    }
  }
  return 2*(2*B - 1 + exp(eps)) / (exp(eps) - 1) * (cnt1 - cnt2/2.);
}

vector<double> HadamardResponse::estimate_all_freqs(const vector<Message> &messages) {
  vector<double> ret(K);
  vector< vector<int> > y(B);
  for (int i = 0; i < B; ++i)
    y[i].resize(b);
  for (Message m : messages)
    y[m.read()/b][m.read()%b]++;
  for (int i = 0; i < B; ++i)
    y[i] = Util::hadamard_transform(y[i]);
  double factor = (2*B - 1 + exp(eps)) / (exp(eps) - 1);
  for (int i = 0; i < K; ++i)
    ret[i] = factor * y[i/(b-1)][(i % (b-1)) + 1];
  return ret;
}

