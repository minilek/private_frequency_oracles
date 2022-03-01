#ifndef SUBSET_SELECTION_H
#define SUBSET_SELECTION_H

#include <boost/random.hpp>
#include "private_frequency_oracle.h"

using namespace std;

class SubsetSelection : public PrivateFrequencyOracle {
  int d;
  double alpha, beta, p;
  boost::uniform_int<> unif_int;

public:
  SubsetSelection(int K_, double eps, bool debug_=false, uint32_t seed=boost::random::mt19937::default_seed);
  Message local_randomizer(int x);
  double estimate_freq(int x, const vector<Message> &messages);
  vector<double> estimate_all_freqs(const vector<Message> &messages);
};

#endif // SUBSET_SELECTION_H
