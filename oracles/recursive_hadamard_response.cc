#include "recursive_hadamard_response.h"
#include "util/util.h"

using namespace std;

RecursiveHadamardResponse::RecursiveHadamardResponse(int K_, double eps, int b_, bool debug_, uint32_t seed)
  : PrivateFrequencyOracle(K_, eps, debug_, seed) {
  k = 1; // log_2 K
  while ((1<<k) < K) 
    ++k;
  K = 1<<k; // round up K to the nearest power of 2
  b = min(k, min(b_, (int)ceil(eps*log2(M_E))));
  block_size = 1<<(b-1);
  B = K / block_size;
  RR = RandomizedResponse(1<<b, eps);
  unif_int = boost::uniform_int<>(0, B-1);
  if (debug) cerr << "K,b=" << K << "," << b << endl;
}

Message RecursiveHadamardResponse::local_randomizer(int x) {
  int r = unif_int(rng);
  // x/B tells the block x is in, and parity(r & t) is rth entry of Hadamard transform of e_t
  Message m = Message(RR.local_randomizer(2*(x/B) + !__builtin_parity(r & x%B)).read()*B + r);
  return m;
}

double RecursiveHadamardResponse::estimate_freq(int x, const vector<Message> &messages) {
  return estimate_all_freqs(messages)[x];
}

vector<double> RecursiveHadamardResponse::estimate_all_freqs(const vector<Message> &messages) {
  double factor = (exp(eps) + (1<<b) - 1) / (exp(eps) - 1);
  vector< vector<int> > Emp(B);
  for (int i = 0; i < B; ++i)
    Emp[i] = vector<int>(block_size);
  
  for (Message msg : messages) {
    int m = msg.read();
    int r = m % B;
    m /= B;
    int bit = m&1;
    if (!bit)
      bit = -1;
    m /= 2;
    Emp[r][m] += bit*factor;
  }
  
  for (int r = 0; r < B; ++r)
    Emp[r] = Util::hadamard_transform(Emp[r]);
  
  vector<int> Y(K);
  for (int r = 0; r < B; ++r)
    for (int j = 0; j < block_size; ++j) 
      Y[j*B + r] += B*Emp[r][j];
  Y = Util::hadamard_transform(Y);
  
  vector<double> ret(K);
  for (int i = 0; i < K; ++i)
    ret[i] = 1.0*Y[i] / K;
  return ret;
}
