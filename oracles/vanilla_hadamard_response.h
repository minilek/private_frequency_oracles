#ifndef VANILLA_HADAMARD_RESPONSE_H
#define VANILLA_HADAMARD_RESPONSE_H

#include "private_frequency_oracle.h"

using namespace std;

class VanillaHadamardResponse : public PrivateFrequencyOracle {
  int k;
  double alpha, beta;
  
public:
  VanillaHadamardResponse(int K_, double eps, bool debug_=false, uint32_t seed=boost::random::mt19937::default_seed); 
  Message local_randomizer(int x);
  double estimate_freq(int x, const vector<Message> &messages);
  vector<double> estimate_all_freqs(const vector<Message> &messages);  
};

#endif // VANILLA_HADAMARD_RESPONSE_H
