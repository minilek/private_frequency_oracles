#include "vanilla_hadamard_response.h"
#include "util/util.h"

using namespace std;

VanillaHadamardResponse::VanillaHadamardResponse(int K_, double eps, bool debug_, uint32_t seed) 
  : PrivateFrequencyOracle(K_, eps, debug_, seed) {
  k = 1;
  while ((1<<k) - 1 < K)
    ++k;
  K = (1<<k) - 1;
  // p = Prob(send orthogonal vector) = exp(eps) / (exp(eps) + 1)
  //
  // unbiased estimator constraints for alpha*[[y in orthog(x)]] * beta:
  // alpha*p + beta = 1 
  // alpha/2 + beta = 0
  //
  // ==> beta = 1 - p*alpha
  // ==> beta = -alpha/2 
  alpha = 2. / (2*exp(eps)/(exp(eps)+1) - 1);
  beta = 1 - exp(eps)/(exp(eps)+1)*alpha;
  if (debug) cerr << "VanillaHadamardResponse K,alpha,beta=" << K << "," << alpha << "," << beta << endl;
}

Message VanillaHadamardResponse::local_randomizer(int x) {
  x++;
  vector<int> bits(k);
  int msb = -1, ret = 0, dotprod_parity = 0;
  for (int i = 0; i < k; ++i) {
    bits[i] = !!(x & (1<<i));
    if (bits[i] == 1)
      msb = i;
  }
  assert(msb >= 0);
  for (int i = 0; i < k; ++i)
    if (i == msb) {
      if(unif_fraction(rng) <= exp(eps) / (exp(eps) + 1)) {
	if (dotprod_parity)
	  ret += 1<<i, dotprod_parity ^= bits[i];
      } else {
	if (!dotprod_parity)
	  ret += 1<<i, dotprod_parity ^= bits[i];
      }
    } else if (unif_fraction(rng) <= 0.5)
      ret += 1<<i, dotprod_parity ^= bits[i];
  return Message(ret);
}

double VanillaHadamardResponse::estimate_freq(int x, const vector<Message> &messages) {
  ++x;
  int cnt = 0;
  for (Message m : messages) {
    int u = m.read(), dotprod_parity = 0;
    for (int i = 0; i < k; ++i)
      if ((x & (1<<i)) && (u & (1<<i)))	  
	dotprod_parity ^= 1;
    if (!dotprod_parity)
      ++cnt;
  }
  return alpha*cnt + beta*messages.size();
}

vector<double> VanillaHadamardResponse::estimate_all_freqs(const vector<Message> &messages) {
  vector<double> ret(K);
  vector<int> y(K + 1), z;
  for (Message m : messages)
    y[m.read()]++;
  z = Util::hadamard_transform(y);
  // if oc = #orthog messages and noc is #non-orthog, then
  // at this point z[i+1] = oc - noc = oc - (n - oc) = 2oc - n
  for (int i = 0; i < K; ++i)
    ret[i] = alpha*((z[i+1] + messages.size())/2.) + beta*messages.size();
  return ret;
}

