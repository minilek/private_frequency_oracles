#ifndef PROJECTIVE_GEOMETRY_H
#define PROJECTIVE_GEOMETRY_H

#include <boost/random.hpp>
#include <vector>
#include "private_frequency_oracle.h"

using namespace std;

class ProjectiveGeometryResponse : public PrivateFrequencyOracle {
  int t, q, cint, cset;
  double p, alpha, beta;
  vector<int> qpows, qinv;
  boost::uniform_int<> unif_int;
  
  int count_orthogonal_messages(const vector<int> &v, int dp_so_far, int at, int index,
				int last_nz, const vector<int> &y);

public:
  ProjectiveGeometryResponse(bool debug_=false, uint32_t seed=boost::random::mt19937::default_seed);
  ProjectiveGeometryResponse (int K_, double eps, int q_, bool debug_=false, uint32_t seed=boost::random::mt19937::default_seed);
  int get_cset();
  int get_cint();
  int get_t();
  Message local_randomizer(int x);  
  vector<int> dp_bottom_up(vector<int> &y); // y[u] is # messages received equal to u
  double estimate_freq(int x, const vector<Message> &messages);
  vector<double> estimate_all_freqs(const vector<Message> &messages);  
};

#endif // PROJECTIVE_GEOMETRY_H
