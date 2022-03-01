#ifndef RECURSIVE_HADAMARD_RESPONSE_H
#define RECURSIVE_HADAMARD_RESPONSE_H

#include <boost/random.hpp>
#include <vector>
#include "private_frequency_oracle.h"
#include "randomized_response.h"

using namespace std;

class RecursiveHadamardResponse : public PrivateFrequencyOracle {
  int b, B, block_size, k;
  RandomizedResponse RR;
  boost::uniform_int<> unif_int;
  
public:
  RecursiveHadamardResponse(int K_, double eps, int b_, bool debug_=false, uint32_t seed=boost::random::mt19937::default_seed);
  Message local_randomizer(int x);
  double estimate_freq(int x, const vector<Message> &messages);
  vector<double> estimate_all_freqs(const vector<Message> &messages);
};

#endif // RECURSIVE_HADAMARD_RESPONSE_H
