#ifndef PRIVATE_FREQUENCY_ORACLE_H
#define PRIVATE_FREQUENCY_ORACLE_H

#include <boost/random.hpp>
#include "util/message.h"
#include "util/util.h"

using namespace std;

class PrivateFrequencyOracle {
  
protected:
  int K;
  double eps;
  bool debug;
  boost::random::mt19937 rng;
  boost::random::uniform_real_distribution<double> unif_fraction;
  
public:
  PrivateFrequencyOracle(bool debug_=false, uint32_t seed=boost::random::mt19937::default_seed);
  PrivateFrequencyOracle(int K_, double eps_, bool debug_, uint32_t seed=boost::random::mt19937::default_seed);
  virtual int universe_size();        
  virtual Message local_randomizer(int x);
  virtual double estimate_freq(int x, const vector<Message> &messages);  
  virtual vector<double> estimate_all_freqs(const vector<Message> &messages); 
};

#endif // PRIVATE_FREQUENCY_ORACLE_H
