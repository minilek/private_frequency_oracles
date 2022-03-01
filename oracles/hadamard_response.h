#ifndef HADAMARD_RESPONSE_H
#define HADAMARD_RESPONSE_H

#include "private_frequency_oracle.h"

using namespace std;

class HadamardResponse : public PrivateFrequencyOracle {
  int b, B, logb;
  double p;
  boost::uniform_int<> unif_int, unif_intb, unif_intB;
  
public:
  HadamardResponse(int K_, double eps, bool debug=false, uint32_t seed=boost::random::mt19937::default_seed); 
  Message local_randomizer(int x);
  double estimate_freq(int x, const vector<Message> &messages);
  vector<double> estimate_all_freqs(const vector<Message> &messages);  
};

#endif // HADAMARD_RESPONSE_H
