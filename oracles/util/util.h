#ifndef PFO_UTIL_H
#define PFO_UTIL_H

using namespace std;

class Util {

  static void generate_divisors(const vector< pair<int,int> > &prime_factorization, vector<int> &ret, int at, int sofar);
  
  // calculate generator for multiplicative group of F_q, q prime
  static int primitive_root(int q, boost::random::mt19937 &rng);

public:
  static vector<int> hadamard_transform(vector<int> x);

  static int fast_mod_pow(int x, int n, int q);
  
  static vector<int> compute_inverse_table(int q, boost::random::mt19937 &rng);
  
  static int first_non_zero(const vector<int> &v);

  // normal a nonzero vector so that it's canonical (first non-entry is a 1)
  static void canonicalize(vector<int> &v, int q, const vector<int> &qinv);

  static vector<int> canonicalized_index_to_vec(int j, int l, int q);
         
  static int canonicalized_vec_to_index(const vector<int> &v, int q, const vector<int> &qpows);

  // vindex is index (1-based indexing) to a canonical vector v in F_q^l
  // returns a 3dim-vector (v[0], g, vsuff_index), where g is the value of the
  // first non-zero entry of vsuffix (vector v with v[0] removed) or 1 if
  // vsuffix is the 0 vector. vsuff_index is the index of vsuffix amongst all
  // canonical vectors in F_q^{l-1}
  static vector<int> decompose_canonical_vector(int vindex, int l, int q,
						const vector<int> &qpows, const vector<int> &qinv);

  static int mod_dot_product(const vector<int> &u, const vector<int> &v, int q);
  
  static bool is_prime(int q);
  
  // pick a random nonzero vector u in F_q^t amongst those with <u,v> = s
  // we assume here that v is not zero (else this isn't possible for all s)
  static vector<int> randvec_with_specified_dot_product(const vector<int>& v, int s, int q,
							const vector<int> &qinv, boost::random::mt19937 &rng);
  
};

#endif // PFO_UTIL_H
