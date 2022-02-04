#include "private_frequency_oracle.h"
#include <vector>

using namespace std;  

PrivateFrequencyOracle::PrivateFrequencyOracle(bool debug_, uint32_t seed) {
  debug = debug_;
  rng = boost::random::mt19937(seed);
}
  
PrivateFrequencyOracle::PrivateFrequencyOracle(int K_, double eps_, bool debug_=false, uint32_t seed) {
    K = K_;
    eps = eps_;
    unif_fraction = boost::random::uniform_real_distribution<double>(0.0, 1.0);
    debug = debug_;
    rng = boost::random::mt19937(seed);
  }

int PrivateFrequencyOracle::universe_size() {
  return K;
}
        
Message PrivateFrequencyOracle::local_randomizer(int x) {
  cerr << "PFO local_randomizer should never be called" << endl;
  cerr.flush();
  assert(false);
  return Message();
}

double PrivateFrequencyOracle::estimate_freq(int x, const vector<Message> &messages) {
  cerr << "PFO estimate_freq should never be called" << endl;
  cerr.flush();
  assert(false);
  return 0.0;
}
  
vector<double> PrivateFrequencyOracle::estimate_all_freqs(const vector<Message> &messages) {
  cerr << "PFO estimate_all_freqs should never be called" << endl;
  cerr.flush();
  assert(false);
  vector<double> ret = vector<double>(K, 0);
  return ret;
}
