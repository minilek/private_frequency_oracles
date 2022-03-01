#include <bits/stdc++.h>
#include <boost/program_options.hpp>
#include <boost/random.hpp>
#include "oracles/util/util.h"
#include "oracles/randomized_response.h"
#include "oracles/rappor.h"
#include "oracles/pi_rappor.h"
#include "oracles/subset_selection.h"
#include "oracles/hadamard_response.h"
#include "oracles/vanilla_hadamard_response.h"
#include "oracles/recursive_hadamard_response.h"
#include "oracles/projective_geometry.h"
#include "oracles/hybrid_projective_geometry.h"

namespace po = boost::program_options;

using namespace std;

int n, trials, K;
int var_rep; // in run_variance_experiments, how many times to run the experiment per value of eps
double var_eps_min, var_eps_max;
int num_var_epsilons;
int kmin, kmax, numk;
int hpgq;
bool debug;
bool rtime, rmse, rvar;
bool run_ss, run_rappor, run_pirappor, run_rr, run_hr, run_vhr, run_rhr, run_pg, run_hpg;
double epsilon;
string base_filename, output_dir;
shared_ptr<istream> input;
ofstream ovar, otime, omse, omax;
vector<int> hist, data;

// if option o1 is set, then at least one of o2 or o3 flags must be true
// also, if either o2 or o3 is true, then o1 must be set
void option_dependency2(const po::variables_map &vm, const char* o1, const char* o2,
			const char* o3, bool o2_triggered, bool o3_triggered) {
  if (o2_triggered && !vm.count(o1))
    throw logic_error(string("Option ' ") + o2 + "' requires option '" + o1 + "'.");
  if (o3_triggered && !vm.count(o1))
    throw logic_error(string("Option ' ") + o3 + "' requires option '" + o1 + "'.");
  if (vm.count(o1) && !(o2_triggered || o3_triggered))
    throw logic_error(string("Option ' ") + o1 + "' requires either option '" + o2 + "' or '" + o3 + "'.");
}

// if option o1 is set, then boolean flag o2 must be true; if o2 is true, o1 must be set
void option_dependency(const po::variables_map &vm, const char* o1, const char* o2, bool o2_triggered) {
  if (o2_triggered && !vm.count(o1))
    throw logic_error(string("Option ' ") + o2 + "' requires option '" + o1 + "'.");
  if (vm.count(o1) && !o2_triggered)
    throw logic_error(string("Option ' ") + o1 + "' requires option '" + o2 + "'.");
}

void init_oracles(PrivateFrequencyOracle **RR, PrivateFrequencyOracle **AsymmRAPPOR,
		  PrivateFrequencyOracle **PIAsymmRAPPOR, PrivateFrequencyOracle **SS,
		  PrivateFrequencyOracle **HR, PrivateFrequencyOracle **VHR,
		  PrivateFrequencyOracle **RHR, PrivateFrequencyOracle **PG,
		  PrivateFrequencyOracle **HPG, vector<PrivateFrequencyOracle*> &algos,
		  vector<string> &algo_names, double eps, int K_, int seed=boost::random::mt19937::default_seed) {
  if (debug) cerr << "initializing oracles" << endl;
  
  if (run_rr) *RR = new RandomizedResponse(K_, eps, debug, seed);
  if (run_rappor) *AsymmRAPPOR = new AsymmetricRAPPOR(K_, eps, debug, seed);
  if (run_pirappor) *PIAsymmRAPPOR = new PIAsymmetricRAPPOR(K_, eps, debug, seed);
  if (run_ss) *SS = new SubsetSelection(K_, eps, debug, seed);
  if (run_hr) *HR = new HadamardResponse(K_, eps, debug, seed);
  if (run_vhr) *VHR = new VanillaHadamardResponse(K_, eps, debug, seed);
  if (run_rhr) *RHR = new RecursiveHadamardResponse(K_, eps, INT_MAX, debug, seed); // b=K will certainly be overwritten
  int q = int(ceil(exp(eps) + 1));
  if (q > 2) {
    if (q % 2 == 0)
      q += 1;
    while (!Util::is_prime(q))
      q += 2;
  }
  if (run_pg) *PG = new ProjectiveGeometryResponse(K_, eps, q, debug, seed); 
  if (run_hpg) *HPG = new HybridProjectiveGeometryResponse(K_, eps, max(2, (int)ceil((exp(eps)+1)/hpgq)), hpgq, debug, seed);

  if (run_rr)
    algos.push_back(*RR), algo_names.push_back("RR");
  if (run_rappor)
    algos.push_back(*AsymmRAPPOR), algo_names.push_back("RAPPOR");
  if (run_pirappor)
    algos.push_back(*PIAsymmRAPPOR), algo_names.push_back("PI-RAPPOR");
  if (run_ss)
    algos.push_back(*SS), algo_names.push_back("SS");
  if (run_hr)
    algos.push_back(*HR), algo_names.push_back("HR");
  if (run_vhr)
    algos.push_back(*VHR), algo_names.push_back("VHR");
  if (run_rhr)
    algos.push_back(*RHR), algo_names.push_back("RHR");
  if (run_pg)
    algos.push_back(*PG), algo_names.push_back("PG");
  if (run_hpg)
    algos.push_back(*HPG), algo_names.push_back(string("HPG") + to_string(hpgq));
  
  if (debug) cerr << "done initializing oracles" << endl;
}

void destroy_oracles(PrivateFrequencyOracle *RR, PrivateFrequencyOracle *AsymmRAPPOR,
		     PrivateFrequencyOracle *PIAsymmRAPPOR, PrivateFrequencyOracle *SS,
		     PrivateFrequencyOracle *HR, PrivateFrequencyOracle *VHR,
		     PrivateFrequencyOracle *RHR, PrivateFrequencyOracle *PG,
		     PrivateFrequencyOracle *HPG) {
  if (run_rr) delete RR;
  if (run_rappor) delete AsymmRAPPOR;
  if (run_pirappor) delete PIAsymmRAPPOR;
  if (run_ss) delete SS;
  if (run_hr) delete HR;
  if (run_vhr) delete VHR;
  if (run_rhr) delete RHR;
  if (run_pg) delete PG;
  if (run_hpg) delete HPG;
}
 
void run_timing_experiments() {
  if (debug) cerr << "==================\nRUNNING TIMING EXPERIMENTS\n==================" << endl;
  PrivateFrequencyOracle *RR, *AsymmRAPPOR, *PIAsymmRAPPOR, *SS, *HR, *VHR, *RHR, *PG, *HPG;
  vector<PrivateFrequencyOracle*> algos;
  vector<string> algo_names;
  int num_algos = run_ss + run_rappor + run_pirappor + run_rr
    + run_hr + run_vhr + run_rhr + run_pg + run_hpg;
  vector< vector<double> > results(num_algos);
  // the actual user data shouldn't affect runtime for any of our oracle implementations
  // so just let each user hold the item "0"
  vector<int> timing_data(n, 0);
  for (int t = 0; t < numk; ++t) {
    int K_ = kmin + (numk == 1 ? 0 : static_cast<int64_t>(t)*(kmax - kmin) / (numk - 1));
    algos.clear();
    algo_names.clear();
    init_oracles(&RR, &AsymmRAPPOR, &PIAsymmRAPPOR, &SS, &HR, &VHR, &RHR, &PG, &HPG, algos, algo_names, epsilon, K_);
    for (int i = 0; i < num_algos; ++i) {
      if (debug) cerr << "K=" << K_ << ", " << algo_names[i] << endl;
      vector<Message> msgs = vector<Message>(n);
      for (int j = 0; j < timing_data.size(); ++j)  
	msgs[j] = algos[i]->local_randomizer(timing_data[j]);
      auto start = chrono::high_resolution_clock::now();
      vector<double> ret = algos[i]->estimate_all_freqs(msgs);
      auto stop = chrono::high_resolution_clock::now();
      if (debug) {
	cerr << "printing first 10 estimates" << endl;
	for (int j = 0; j < 10; ++j)
	  cerr << ret[j] << " ";
	cerr << endl;
      }
      auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
      results[i].push_back(duration.count());
    }
    destroy_oracles(RR, AsymmRAPPOR, PIAsymmRAPPOR, SS, HR, VHR, RHR, PG, HPG);
  }
  for (int i = 0; i < algos.size(); ++i)
    for (int t = 0; t < numk; ++t) 
      otime << algo_names[i] << "," << kmin + (numk == 1 ? 0 : static_cast<int64_t>(t)*(kmax - kmin) / (numk - 1)) << "," << results[i][t] << endl;
}

void run_variance_experiments() {
  if (debug) cerr << "==================\nRUNNING PRIVACY/UTILITY TRADEOFF EXPERIMENTS\n==================" << endl;
  PrivateFrequencyOracle *RR, *AsymmRAPPOR, *PIAsymmRAPPOR, *SS, *HR, *VHR, *RHR, *PG, *HPG;
  vector<PrivateFrequencyOracle*> algos;
  vector<string> algo_names;
  init_oracles(&RR, &AsymmRAPPOR, &PIAsymmRAPPOR, &SS, &HR, &VHR, &RHR, &PG, &HPG, algos, algo_names, 5.0, 22000); // ugly hack to get algo_names set
  destroy_oracles(RR, AsymmRAPPOR, PIAsymmRAPPOR, SS, HR, VHR, RHR, PG, HPG);
  int num_algos = run_ss + run_rappor + run_pirappor + run_rr
    + run_hr + run_vhr + run_rhr + run_pg + run_hpg;
  vector< vector<double> > results(num_algos);
  double increment = (var_eps_max - var_eps_min) / (num_var_epsilons - 1);
  for (int t = 0; t < num_var_epsilons; ++t) {
    double eps = var_eps_min + t*increment;
    vector<double> mse_averages(num_algos, 0.0);
    if (debug) cerr << "epsilon=" << eps << endl;
    for (int rr = 0; rr < var_rep; ++rr) {
      algos = vector<PrivateFrequencyOracle*>();
      algo_names = vector<string>();
      init_oracles(&RR, &AsymmRAPPOR, &PIAsymmRAPPOR, &SS, &HR, &VHR, &RHR, &PG, &HPG, algos, algo_names, eps, K, boost::random::mt19937::default_seed + rr);
      for (int i = 0; i < num_algos; ++i) {
	if (debug) cerr << algo_names[i] << ":" << endl;
	vector<Message> msgs = vector<Message>(n);
	for (int j = 0; j < data.size(); ++j)  
	  msgs[j] = algos[i]->local_randomizer(data[j]);
	vector<double> ret = algos[i]->estimate_all_freqs(msgs);
	if (debug) {
	  cerr << "printing first 10 estimates (actuals)" << endl;
	  for (int j = 0; j < 10; ++j)
	    cerr << ret[j] << "(" << hist[j] << ") ";
	  cerr << endl;
	}
	for (int j = 0; j < K; ++j) 
	  mse_averages[i] += (ret[j]-hist[j])*(ret[j]-hist[j]);
	if (debug) cerr << "done computing mse averages" << endl;
      }
      destroy_oracles(RR, AsymmRAPPOR, PIAsymmRAPPOR, SS, HR, VHR, RHR, PG, HPG);
    }
    for (int i = 0; i < num_algos; ++i)
      results[i].push_back(mse_averages[i] / (var_rep * K));
  }
  for (int i = 0; i < algos.size(); ++i)
    for (int t = 0; t < num_var_epsilons; ++t)
      ovar << algo_names[i] << "," << var_eps_min + t*increment << "," << results[i][t] << endl;
}
 
void run_mse_experiments() {
  if (debug) cerr << "==================\nRUNNING MSE EXPERIMENTS\n==================" << endl;
  PrivateFrequencyOracle *RR, *AsymmRAPPOR, *PIAsymmRAPPOR, *SS, *HR, *VHR, *RHR, *PG, *HPG;
  vector<PrivateFrequencyOracle*> algos;
  vector<string> algo_names;
  init_oracles(&RR, &AsymmRAPPOR, &PIAsymmRAPPOR, &SS, &HR, &VHR, &RHR, &PG, &HPG, algos, algo_names, epsilon, K);
  
  vector< vector<double> > results(algos.size()), max_results(algos.size());
  
  for (int rr = 0; rr < trials; ++rr)
    for (int i = 0; i < algos.size(); ++i) {
      if (debug) cerr << "trial# " << rr << ", " << algo_names[i] << endl;
      vector<Message> msgs = vector<Message>(n);
      for (int j = 0; j < data.size(); ++j)  
	msgs[j] = algos[i]->local_randomizer(data[j]);
      if (debug) cerr << "estimating all frequencies" << endl;
      vector<double> ret = algos[i]->estimate_all_freqs(msgs);
      if (debug) {
	cerr << "printing first 10 estimates (actuals)" << endl;
	for (int j = 0; j < 10; ++j)
	  cerr << ret[j] << "(" << hist[j] << ") ";
	cerr << endl;
      }
      double bias = 0;
      for (int j = 0; j < K; ++j)
	bias += (ret[j] - hist[j]);
      if (debug) cerr << "avg bias in universe is " << bias/K << endl;
      double mse = 0, maxdev = 0;
      for (int j = 0; j < K; ++j) {
	mse += (ret[j]-hist[j])*(ret[j]-hist[j]);
	maxdev = max(maxdev, abs(ret[j]-hist[j]));
      }
      results[i].push_back(mse / K);
      max_results[i].push_back(maxdev);
    }
  
  for (int i = 0; i < algos.size(); ++i) {
    sort(results[i].begin(), results[i].end());
    sort(max_results[i].begin(), max_results[i].end());
    assert(results[i].size() == max_results[i].size());
    for (int j = 0; j < results[i].size(); ++j) {
      omse << algo_names[i] << "," << results[i][j] << endl;
      omax << algo_names[i] << "," << max_results[i][j] << endl;
    }
  }

  destroy_oracles(RR, AsymmRAPPOR, PIAsymmRAPPOR, SS, HR, VHR, RHR, PG, HPG);
}

void set_base_filename(string s) {
  if (!s.compare("-"))
    base_filename = "pfotest";
  int i = s.find_last_of('/');
  if (i != string::npos)
    s = s.substr(i + 1);
  i = s.find_last_of('.');
  if (i == string::npos)
    base_filename = string(s);
  else
    base_filename = s.substr(0, i);    
}

void read_histogram_data() {
  hist.resize(K);
  data.resize(n);
  fill(hist.begin(), hist.end(), 0);
  int x, cnt = 0;
  while (*input >> x) {
    assert(cnt < n);
    assert(x>=0 && x<K);
    data[cnt++] = x;
    hist[x]++;
  }
  assert(cnt == n);
}

int main(int argc, char* argv[]) {
  try {

    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h",                                                                              "print usage message")
      ("num-users,n",                   po::value<int>(&n)->required(),                       "number of users")
      ("universe-size,u",               po::value<int>(&K),                                   "universe size")
      ("num-trials,t",                  po::value<int>(&trials),                              "number of random trials to run in mse experiments")
      ("epsilon,e",                     po::value<double>(&epsilon),                          "privacy parameter epsilon for mse/timing experiments")
      ("input-file,i",                  po::value<string>()->default_value("-"),              "file path specifying users' items")
      ("output-dir,o",                  po::value<string>()->default_value("pfotmp/"),        "dir in which to write experiment output files")
      ("debug",                         po::bool_switch(&debug)->default_value(false),        "debug mode prints more to cerr")
      ("include-ss,ss",                 po::bool_switch(&run_ss)->default_value(false),       "include SubsetSelection in experiments")
      ("include-rappor,rappor",         po::bool_switch(&run_rappor)->default_value(false),   "include RAPPOR in experiments")
      ("include-pirappor,pirappor",     po::bool_switch(&run_pirappor)->default_value(false), "include PI-RAPPOR in experiments")
      ("include-rr,rr",                 po::bool_switch(&run_rr)->default_value(false),       "include RandomizedResponse in experiments")
      ("include-hr,hr",                 po::bool_switch(&run_hr)->default_value(false),       "include HadamardResponse in experiments")
      ("include-vhr,vhr",               po::bool_switch(&run_vhr)->default_value(false),      "include VanillaHadamardResponse in experiments")
      ("include-rhr,rhr",               po::bool_switch(&run_rhr)->default_value(false),      "include RecursiveHadamardResponse in experiments")
      ("hpgq",                          po::value<int>(&hpgq),                                "value of q to use with HybridProjectiveGeometry")
      ("include-pg,pg",                 po::bool_switch(&run_pg)->default_value(false),       "include ProjectiveGeometry in experiments")
      ("include-hpg,hpg",               po::bool_switch(&run_hpg)->default_value(false),      "include HybridProjectiveGeometry in experiments")
      ("run-variance-experiments,rvar", po::bool_switch(&rvar)->default_value(false),         "run experiments to test variance")
      ("run-timing-experiments,rtime",  po::bool_switch(&rtime)->default_value(false),        "run experiments to test runtimes")
      ("run-mse-experiments,rmse",      po::bool_switch(&rmse)->default_value(false),         "run experiments to test mse")
      ("variance-repetitions,vr",       po::value<int>(&var_rep),                             "number of mse's to average per eps when running variance experiments")
      ("var-epsilon-min,vemn",          po::value<double>(&var_eps_min),                      "smallest eps to use in variance experiment")
      ("var-epsilon-max,vemx",          po::value<double>(&var_eps_max),                      "largest eps to use in variance experiment")
      ("num-var-epsilons,nve",          po::value<int>(&num_var_epsilons),                    "number of epsilons to use in variance experiment")
      ("universe-min,umin",             po::value<int>(&kmin),                                "smallest universe size to use in timing experiments")
      ("universe-max,umax",             po::value<int>(&kmax),                                "largest universe size to use in timing experiments")
      ("num-universes,nu",              po::value<int>(&numk),                                "number of universe sizes to try between umin and umax")
      ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {  
      cout << desc << endl;
      return 0;
    }
    
    po::notify(vm);   

    option_dependency2(vm, "epsilon", "run-mse-experiments", "run-timing-experiments", rmse, rtime);
    option_dependency2(vm, "universe-size", "run-mse-experiments", "run-variance-experiments", rmse, rvar);
    option_dependency(vm, "hpgq", "include-hpg", run_hpg);
    option_dependency(vm, "variance-repetitions", "run-variance-experiments", rvar);
    option_dependency(vm, "var-epsilon-min", "run-variance-experiments", rvar);
    option_dependency(vm, "var-epsilon-max", "run-variance-experiments", rvar);
    option_dependency(vm, "num-var-epsilons", "run-variance-experiments", rvar);
    option_dependency(vm, "universe-min", "run-timing-experiments", rtime);
    option_dependency(vm, "universe-max", "run-timing-experiments", rtime);
    option_dependency(vm, "num-universes", "run-timing-experiments", rtime);
    option_dependency(vm, "num-trials", "run-mse-experiments", rmse);

    if (rvar)
      assert((num_var_epsilons > 1) && (var_eps_min != var_eps_max));
    if (rtime)
      assert(kmin <= kmax && numk >= 1);

    // set input and output files
    if (!vm["input-file"].as<string>().compare("-"))
      input.reset(&cin, [](...){});
    else
      input.reset(new ifstream(vm["input-file"].as<string>()));
    set_base_filename(vm["input-file"].as<string>());

    output_dir = vm["output-dir"].as<string>();
    if (output_dir[output_dir.size() - 1] != '/')
      output_dir += string("/");
    if (rvar) ovar = ofstream(output_dir + base_filename + string(".var"));
    if (rtime) otime = ofstream(output_dir + base_filename + string(".time"));
    if (rmse) {
      omse = ofstream(output_dir + base_filename + string(".mse"));
      omax = ofstream(output_dir + base_filename + string(".max"));
    }

    // make sure at least one oracle is being tested
    if (!(run_ss || run_rappor || run_pirappor || run_rr || run_hr || run_vhr || run_rhr || run_pg || run_hpg))
      throw logic_error("At least one algorithm must be included in experiment.");

    // make sure at least one type of experiment is being run
    if (!(rvar || rtime || rmse))
      throw logic_error("At least one type of experiment must be run.");

    // read in user data for the experiments if not a timing experiment
    if (rmse || rvar)
      read_histogram_data();

    if (rvar) run_variance_experiments();
    if (rtime) run_timing_experiments();
    if (rmse) run_mse_experiments();

  } catch(std::exception& e) {
    cout << e.what() << endl;
    return 1;
  }
  
  return 0;
}
