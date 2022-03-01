#ifndef RAPPOR_H
#define RAPPOR_H

#include "private_frequency_oracle.h"
#include <vector>

using namespace std;

class SymmetricRAPPOR : public PrivateFrequencyOracle {  
  double p;
  
public:
  SymmetricRAPPOR(int K, double eps, bool debug_=false, uint32_t seed=boost::random::mt19937::default_seed);
  Message local_randomizer(int x);
  double estimate_freq(int x, const vector<Message> &messages);
  vector<double> estimate_all_freqs(const vector<Message> &messages);
};

class AsymmetricRAPPOR : public PrivateFrequencyOracle {  
  double p1, p2, alpha, beta;
  
public:
  AsymmetricRAPPOR(int K, double eps, bool debug_=false, uint32_t seed=boost::random::mt19937::default_seed);
  Message local_randomizer(int x);
  double estimate_freq(int x, const vector<Message> &messages);
  vector<double> estimate_all_freqs(const vector<Message> &messages);
};

#endif // RAPPOR_H
