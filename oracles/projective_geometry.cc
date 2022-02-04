#include "projective_geometry.h"
#include "util/util.h"

using namespace std;

int ProjectiveGeometryResponse::count_orthogonal_messages(const vector<int> &v, int dp_so_far, int at, int index,
							int last_nz, const vector<int> &y) {
  if (at == v.size()) {
    if (dp_so_far)
      return 0;
    else
      return y[index];
  } else if (at != last_nz) {
    int ret = 0;
    for (int d = 0; d < q; ++d)
      ret += count_orthogonal_messages(v, (dp_so_far + (int64_t)d*v[at]) % q,
				       at+1, index+d*qpows[t-1-at], last_nz, y);
    return ret;
  } else {
    // need dp_so_far + d*v[last_nz] = 0 mod q
    int d = (((int64_t)q - dp_so_far) * qinv[v[last_nz]]) % q;
    return count_orthogonal_messages(v, (dp_so_far + (int64_t)d*v[at]) % q,
				     at+1, index+d*qpows[t-1-at], last_nz, y);
  }
}

ProjectiveGeometryResponse::ProjectiveGeometryResponse(bool debug_, uint32_t seed_)
  : PrivateFrequencyOracle(debug_, seed_) { }
  
ProjectiveGeometryResponse::ProjectiveGeometryResponse (int K_, double eps, int q_, bool debug_, uint32_t seed_)
  : PrivateFrequencyOracle(K_, eps, debug_, seed_) {
  t = 1;
  q = q_;
  int qpow = q_;
  qpows.push_back(1);
  qpows.push_back(q);
  while ((t < 2) || ((qpow-1) / (q-1) < K)) {
    ++t;
    qpow *= q;
    qpows.push_back(qpow);
  }
  K = (qpow - 1) / (q - 1);
  qinv = Util::compute_inverse_table(q, rng);
  /*  if (debug) {cerr << "printing inverse table" << endl;
    for (int i = 0; i < q; ++i) cerr << qinv[i] << " ";
    cerr << endl;
    }*/
  cint = (qpows[t-2] - 1) / (q - 1);
  cset = (qpows[t-1] - 1) / (q - 1);
  p = 1.0 / ((exp(eps) - 1)*cset + K);
  // alpha*s + beta = 1; alpha*r + beta = 0
  // solves to alpha = 1/(s-r); beta = -alpha*r
  // here s = exp(eps)*p*cset; r = p*((exp(eps)-1)*cint + cset
  double s = exp(eps)*p*cset;
  double r = p*((exp(eps)-1)*cint + cset);
  alpha = 1.0 / (s - r);
  beta = -alpha * r;
  if (debug) cerr << "ProjectiveGeometry K,q,t,cint,cset,alpha,beta=" << K << "," << q << "," << t << "," << cint << "," << cset << "," << alpha << "," << beta << endl;
  unif_int = boost::uniform_int<>(1, q-1); 
}

int ProjectiveGeometryResponse::get_cset() {
  return cset;
}

int ProjectiveGeometryResponse::get_cint() {
  return cint;
}

int ProjectiveGeometryResponse::get_t() {
  return t;
}

Message ProjectiveGeometryResponse::local_randomizer(int x) {
  vector<int> v = Util::canonicalized_index_to_vec(x, t, q), ret;
  if (unif_fraction(rng) <= p*exp(eps)*cset) // pick a random element in S_v (viewed as an element of (F_q)^t)
    ret = Util::randvec_with_specified_dot_product(v, 0, q, qinv, rng);
  else // pick a random element not in S_v (viewed as an element of (F_q)^t)
    ret = Util::randvec_with_specified_dot_product(v, unif_int(rng), q, qinv, rng);
  Util::canonicalize(ret, q, qinv);
  return Message(Util::canonicalized_vec_to_index(ret, q, qpows));
}

// y[u] is # messages received equal to u
// returns length-K vector T, where T[u] = sum_{v in S_u} y_v
// via space-optimized bottom-up dynamic programming
vector<int> ProjectiveGeometryResponse::dp_bottom_up(vector<int> &y) {
  // a is a prefix of u, b is a suffix of v in {0,...,K-1} (but thought of as
  // a canonical vec). DP: f(a, b, z) = sum_{u starting with a : <suffix(u),b> = z} y_u
  // at the end we care about f(empty, v, 0) (base case is when b=empty)
  // table entries with a having length l only depend on table entries
  // with a having length l+1, so can save memory by doing bottom-up DP
  // and reusing memory. last is the dp array for length+1, and next is for length
  // ultimately we want the answer for length==0, and base
  // case is length==t
  
  // number of possible (a,b) combinations for any fixed value of length=l
  // is K+1 = (q^t-1)/(q-1) + 1 for l==t, where the "+ 1" is for a == 0. Then 
  // for smaller length l it is ((q^l-1)/(q-1) + 1) x (q^{t-l}) = (#a's) x (#b's)
  // Now, put answers in last for the base cases (l==t)
  // N is the biggest size DP table we need for any fixed value of l
  assert(y.size() == K);
  int N = K + 1; // mem table requirement for l==0, since we don't store answers for z!=0 to save memory
  for (int l = 1; l < t; ++l) 
    N = max(N, ((qpows[l]-1)/(q-1) + 1) * ((qpows[t-l]-1)/(q-1) + 1) * q);
  vector<int> last(N), next(N);
  
  // answers for base cases when l==t
  // answers for a==0 should be 0; only nonzero answers are when
  // a!=0 and z==0; note non-zero a must be canonical when l==t
  for (int a = 1; a <= K; ++a)
    last[a] = y[a-1];   
  
  int lastA = K+1, lastB = 1, curA = 0, curB = 0;
  vector<int> ret(K);
  
  // get answers for non base cases
  for (int length = t - 1; length >= 0; --length) {
    curA = (qpows[length] - 1) / (q-1) + 1, curB = (qpows[t - length] - 1) / (q-1) + 1;
    fill(next.begin(), next.end(), 0);
    for (int b = 0; b < curB; ++b) {
      vector<int> decomp = Util::decompose_canonical_vector(b, t - length, q, qpows, qinv);
      int vb0 = decomp[0], ginv = qinv[decomp[1]], vbsuff_index = decomp[2];
      for (int a = 0; a < curA; ++a) {
	if (!length) { // only care about z==0
	  int calc = last[vbsuff_index*lastA*q + 0*q + 0]; // append 0 to a
	  calc += last[vbsuff_index*lastA*q + 1*q + (((int64_t)q - vb0) * ginv) % q]; // append 1 to a
	  next[b] = calc;
	} else {
	  int extension = a ? (2 + (a-1)*q) : 0; // index of a0 in last ('2 +' for a=0 and [0,..,0,1], then each previous nonzero a has q extensions, unless a==0)
	  for (int z = 0; z < q; ++z) {
	    int calc = 0;
	    for (int d = 0; d <= (a ? q-1 : 1); ++d) { // extend a by digit d; cannot choose d = 2,3,..,q-1 if a==0
	      int new_dot_prod = ((((int64_t)q + z - vb0*d) % q) * ginv) % q; // what new dot prod needs to be after extending a by d
	      if (length == t-1)  // to save memory we stored last differently in this case
		calc += (new_dot_prod ? 0 : last[extension + d]);
	      else 
		calc += last[vbsuff_index*lastA*q + (extension+d)*q + new_dot_prod];
	    }
	    next[b*curA*q + a*q + z] = calc;
	  }
	}
      }
    }
    swap(last, next);
    lastA = curA;
    lastB = curB;
  }
  for (int i = 0; i < K; ++i)
    ret[i] = last[i + 1];
  return ret;
}

double ProjectiveGeometryResponse::estimate_freq(int x, const vector<Message> &messages) {
  vector<int> v = Util::canonicalized_index_to_vec(x, t, q), u;
  int cnt = 0; // count number of messages orthogonal to v
  for (Message m : messages) {
    u = Util::canonicalized_index_to_vec(m.read(), t, q);
    if (Util::mod_dot_product(u, v, q) == 0)
      cnt++;
  }
  return alpha*cnt + beta*messages.size();
}

vector<double> ProjectiveGeometryResponse::estimate_all_freqs(const vector<Message> &messages) {
  vector<double> ret(K);
  vector<int> y(K);
  for (Message m : messages) 
    y[m.read()]++;
  if (t <= 3) {
    // O(k^2 / q) time; faster when K <= q^2 t; since K = (q^t-1)/(q-1) so that q ~ K^{1/(t-1)},
    // this happens roughly when K <= t K^{2/(t-1)}, which means t <= 3 (so either t=2 or t=3)
    for (int x = 0; x < K; ++x) {
      int last_nz = -1;
      vector<int> v = Util::canonicalized_index_to_vec(x, t, q);
      for (int i = v.size() - 1; i >= 0; --i)
	if (v[i]) {
	  last_nz = i;
	  break;
	}
      assert(last_nz != -1);
      int offset = 0;
      for (int j = t - 1; j >= 0; --j)  {
	// sum y_u for orthogonal u whose first nonzero is in coord j
	ret[x] += count_orthogonal_messages(v, v[j], j + 1, offset, last_nz, y);
	offset += qpows[t - 1 - j]; // number of canonical vectors whose first nonzero is in coord j
      }
      ret[x] = alpha*ret[x] + beta*messages.size();
    }
  } else {
    // O(kqt) time
    vector<int> T = dp_bottom_up(y); // results from dynamic programming
    for (int i = 0; i < K; ++i)
      ret[i] = alpha*T[i] + beta*messages.size();
  }
  return ret;
}
