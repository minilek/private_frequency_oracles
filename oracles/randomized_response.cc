#include "randomized_response.h"
#include <vector>
using namespace std;

RandomizedResponse::RandomizedResponse(bool debug_, uint32_t seed_)
  : PrivateFrequencyOracle(debug_, seed_) { }

RandomizedResponse::RandomizedResponse(int K, double eps, bool debug_, uint32_t seed_)
  : PrivateFrequencyOracle(K, eps, debug_, seed_) {
  p = exp(eps) / (exp(eps) + K - 1);
  unif_int = boost::uniform_int<>(0, K-2);
}

Message RandomizedResponse::local_randomizer(int x) {
  if (unif_fraction(rng) <= p)
    return Message(x);
  else {
    int r = unif_int(rng); 
    if (r < x)
      return Message(r);
    else 
      return Message(r+1);
  }
}
            
double RandomizedResponse::estimate_freq(int x, const vector<Message> &messages) {
  int cnt = 0;
  for (Message m : messages)
    if (m.read() == x)
      cnt++;
  double alpha = (exp(eps)+ K - 1) / (exp(eps) - 1);
  double beta = -1.0 / (exp(eps) - 1);
  return alpha*cnt + beta*messages.size();
}

vector<double> RandomizedResponse::estimate_all_freqs(const vector<Message> &messages) {
  vector<double> cnts = vector<double>(K, 0.0);
  for (Message m : messages) 
    cnts[m.read()]++;    
  double alpha = (exp(eps)+ K - 1) / (exp(eps) - 1);
  double beta = -1.0 / (exp(eps) - 1);
  for (int i = 0; i < K; ++i)
    cnts[i] = alpha*cnts[i] + beta*messages.size();
  return cnts;
}
