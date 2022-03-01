#include "rappor.h"

using namespace std;

SymmetricRAPPOR::SymmetricRAPPOR(int K, double eps, bool debug_, uint32_t seed)
  : PrivateFrequencyOracle(K, eps, debug_, seed) {
  p = 1.0 / (exp(eps) + 1);
}

Message SymmetricRAPPOR::local_randomizer(int x) {
  vector<int> ret(K, 0);
  ret[x] = 1;
  for (int i = 0; i < K; ++i) 
    if (unif_fraction(rng) <= p)
      ret[i] = 1 - ret[i];
  return Message(ret);
}

double SymmetricRAPPOR::estimate_freq(int x, const vector<Message> &messages) {
  int cnt = 0;
  for (Message m : messages)
    if (m[x] == 1)
      ++cnt;
  double alpha = 1.0 / (1.0 - 2*p);
  double beta = -p / (1.0 - 2*p);
  return alpha*cnt + beta*messages.size();
}

vector<double> SymmetricRAPPOR::estimate_all_freqs(const vector<Message> &messages) {
  vector<double> cnts(K, 0);
  for (int x = 0; x < K; ++x)
    for (Message m : messages)
      if (m[x] == 1)
	cnts[x]++;
  double alpha = 1.0 / (1.0 - 2*p);
  double beta = -p / (1.0 - 2*p);
  for (int x = 0; x < K; ++x)
    cnts[x] = alpha*cnts[x] + beta*messages.size();
  return cnts;
}

AsymmetricRAPPOR::AsymmetricRAPPOR(int K, double eps, bool debug_, uint32_t seed)
  : PrivateFrequencyOracle(K, eps, debug_, seed) {
  p1 = 1.0 / (exp(eps) + 1);
  p2 = 0.5;
  alpha = 1.0 / (1.0 - p1 - p2);
  beta = -p1 / (1.0 - p1 - p2);
  if (debug) cerr << "AssymetricRAPPOR alpha,beta=" << alpha << "," << beta << endl;
}

Message AsymmetricRAPPOR::local_randomizer(int x) {
  vector<int> ret(K, 0);
  ret[x] = 1;
  for (int i = 0; i < K; ++i)
    if (i != x) {
      if (unif_fraction(rng) <= p1)
	ret[i] = 1;
    } else {
      if (unif_fraction(rng) <= p2)
	ret[i] = 0;
    }
  return Message(ret);
}

double AsymmetricRAPPOR::estimate_freq(int x, const vector<Message> &messages) {
  int cnt = 0;
  for (Message m : messages)
    if (m[x] == 1)
      ++cnt;
  return alpha*cnt + beta*messages.size();
}

vector<double> AsymmetricRAPPOR::estimate_all_freqs(const vector<Message> &messages) {
  vector<double> cnts(K, 0);
  for (int x = 0; x < K; ++x)
    for (Message m : messages)
      if (m[x] == 1)
	cnts[x]++;
  for (int x = 0; x < K; ++x)
    cnts[x] = alpha*cnts[x] + beta*messages.size();
  return cnts;
}
