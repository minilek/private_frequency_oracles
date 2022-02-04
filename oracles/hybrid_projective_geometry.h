#ifndef HYBRID_PROJECTIVE_GEOMETRY_H
#define HYBRID_PROJECTIVE_GEOMETRY_H

#include "projective_geometry.h"

class HybridProjectiveGeometryResponse : public PrivateFrequencyOracle {
  int h, q, t, b, cset, cint;
  double p, alpha, beta, gamma;
  ProjectiveGeometryResponse PG;
  vector<int> qinv, qpows;
  boost::uniform_int<> unif_int, unif_inth, unif_intq;

public:
  HybridProjectiveGeometryResponse(int K_, double eps_, int h_, int q_, bool debug_=false, uint32_t seed=boost::random::mt19937::default_seed);
  Message local_randomizer(int x);
  double estimate_freq(int x, const vector<Message> &messages);
  vector<double> estimate_all_freqs(const vector<Message> &messages);
};

#endif // HYBRID_PROJECTIVE_GEOMETRY_H
