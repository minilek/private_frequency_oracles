#ifndef PI_RAPPOR_H
#define PI_RAPPOR_H

#include <boost/random.hpp>
#include <vector>
#include "private_frequency_oracle.h"

using namespace std;

class PIAsymmetricRAPPOR : public PrivateFrequencyOracle {  
  double p1, p2, alpha, beta;
  int q, t, zero_threshold;
  vector<int> qpows, qinv;
  boost::uniform_int<> unif_intq, unif_int0, unif_int1;

  vector<int> num_to_vec(int x, int l);
  int vec_to_num(const vector<int> &v);
  
public:
  PIAsymmetricRAPPOR(int K_, double eps, bool debug=false, uint32_t seed=boost::random::mt19937::default_seed);
  Message local_randomizer(int x);
  double estimate_freq(int x, const vector<Message> &messages);
  vector<int> dp_bottom_up(const vector<int> &y); // y[u] is the number of messages received equal to u
  vector<double> estimate_all_freqs(const vector<Message> &messages);
};

#endif // PI_RAPPOR_H
