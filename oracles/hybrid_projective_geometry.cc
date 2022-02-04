#include "hybrid_projective_geometry.h"
#include "util/util.h"

using namespace std;

/*  
  For each u in S_v, PG sends it with probability exp(eps)*p where p = 1 / ((exp(eps) - 1)*cset + b)
  -- what's in the PDF: some items have probability exp(eps)*p of being output, others have probability p
  -- i can view my input x as a pair (B, v) (B is a block ID in {0,..,h-1}, and v is in {0,...,b-1} for b = k/h)
     want to send (block, random u s.t. <u,v> = 0)
  -- Randomized response: sends (B, PG output) w.p. exp(eps)*p'
                          (this means (B, u) for any fixed u in S_v     is sent w.p. exp(eps)*p)
                          (and        (B, u') for any fixed u notin S_v is sent w.p. 
                          else sends (B', z) for random b'!=b and unif random z
 */
HybridProjectiveGeometryResponse::HybridProjectiveGeometryResponse(int K_, double eps_, int h_, int q_, bool debug_, uint32_t seed_)
  : PrivateFrequencyOracle(K_, eps_, debug_, seed_) {
  h = h_;
  q = q_;
  t = 0;
  int qpow = 1;
  qpows.push_back(qpow);
  while ((t < 2) || ((qpow - 1) / (q - 1) * h < K)) {
    qpow *= q;      
    ++t;
    qpows.push_back(qpow);
  }
  b = (qpow - 1) / (q - 1);
  PG = ProjectiveGeometryResponse(b, eps, q);
  assert(PG.universe_size() == b);
  assert(PG.get_t() == t);
  cset = PG.get_cset();
  cint = PG.get_cint();
  qinv = Util::compute_inverse_table(q, rng);
  K = b*h;
  if (debug) cerr << "h,q,t,b,K " << h << "," << q << "," << t << "," << b << "," << K << endl;
  unif_int = boost::uniform_int<>(0, PG.universe_size() - 1);
  unif_inth = boost::uniform_int<>(0, h-2);
  unif_intq = boost::uniform_int<>(1, q-1);
  p = 1.0 / (b*h + (exp(eps) - 1)*cset);
  alpha = (1.0*b*h + (exp(eps)-1)*cset) / ((exp(eps)-1)*(cset - cint));
  beta = -(1.0*b*h + (exp(eps)-1)*cset) / ((exp(eps)-1)*(cset - cint)) * cint / cset;
  gamma = -alpha*p*cset - beta*p*b;
  if (debug) cerr << "alpha,beta,gamma " << alpha << "," << beta << "," << gamma << endl;
}

Message HybridProjectiveGeometryResponse::local_randomizer(int x) {
  int block, vec;
  if (unif_fraction(rng) <= exp(eps)*p*cset) {
    block = x / b;
    vector<int> v = Util::canonicalized_index_to_vec(x % b, t, q);
    vector<int> v2 = Util::randvec_with_specified_dot_product(v, 0, q, qinv, rng);
    Util::canonicalize(v2, q, qinv);
    vec = Util::canonicalized_vec_to_index(v2, q, qpows);
  } else {
    int A = (h-1) * b; // #universe elts with different block
    int B = b - cset; // #universe elts with same block, but 2nd coord not in S_v
    if (unif_fraction(rng) <= 1. * A / (A + B)) {
      block = unif_inth(rng);
      if (block >= x / b)
	++block;
      vec = unif_int(rng);
    } else {
      block = x / b;
      vector<int> v = Util::canonicalized_index_to_vec(x % b, t, q);
      vector<int> v2 = Util::randvec_with_specified_dot_product(v, unif_intq(rng), q, qinv, rng);
      Util::canonicalize(v2, q, qinv);
      vec = Util::canonicalized_vec_to_index(v2, q, qpows);
    }
    
  }
  return Message(block*b + vec);
}

double HybridProjectiveGeometryResponse::estimate_freq(int x, const vector<Message> &messages) {
  vector<int> v = Util::canonicalized_index_to_vec(x % b, t, q), u;
  int block = x / b, bl;
  int cnt1 = 0, cnt2 = 0; // the first two sums in the PDF
  for (Message m : messages) {
    bl = m.read() / b;
    u = Util::canonicalized_index_to_vec(m.read() % b, t, q);
    int mdp = Util::mod_dot_product(u, v, q);
    if (bl==block) {
      ++cnt2;
      if (!mdp)
	++cnt1;
    }
  }
  return alpha*cnt1 + beta*cnt2 + gamma*messages.size();
}

vector<double> HybridProjectiveGeometryResponse::estimate_all_freqs(const vector<Message> &messages) {
  vector< vector<int> > cnts(h);
  vector<int> cnts2(h);
  for (int i = 0; i < h; ++i)
    cnts[i] = vector<int>(b);
  for (Message m : messages) {
    int block = m.read() / b, val = m.read() % b;
    assert(block < h);
    cnts[block][val]++;
    cnts2[block]++;
  }
  vector< vector<int> > rets(h);
  for (int block = 0; block < h; ++block)
    rets[block] = PG.dp_bottom_up(cnts[block]);
  vector<double> ret(K);
  for (int i = 0; i < K; ++i)
    ret[i] = alpha*rets[i / b][i % b] + beta*cnts2[i / b] + gamma*messages.size();
  return ret;
}
