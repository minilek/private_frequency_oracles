#ifndef RANDOMIZED_RESPONSE_H
#define RANDOMIZED_RESPONSE_H

#include "private_frequency_oracle.h"
#include <boost/random.hpp>

class RandomizedResponse : public PrivateFrequencyOracle {  
  double p;
  boost::uniform_int<> unif_int;
  
public:
  RandomizedResponse(bool debug_=false, uint32_t seed=boost::random::mt19937::default_seed);
  RandomizedResponse(int K, double eps, bool debug_=false, uint32_t seed_=boost::random::mt19937::default_seed);
  Message local_randomizer(int x);
  double estimate_freq(int x, const vector<Message> &messages);  
  vector<double> estimate_all_freqs(const vector<Message> &messages); 
};

#endif // RANDOMIZED_RESPONSE_H
