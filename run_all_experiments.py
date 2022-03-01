#!/usr/bin/env python3

import argparse
import bisect
import datetime
import os
from pathlib import Path
import random
import shutil
import sys

from absl import app, flags

from pfoparams import Params

FLAGS = flags.FLAGS

flags.DEFINE_string("output_dir", "pfotmp", "output directory experimental results")
flags.DEFINE_boolean("debug", False, "produces debugging output")
    
def print_data(filename, data):    
    filepath = dirpath / 'data' / filename
    with filepath.open('w', encoding='utf-8') as f:
        print('printing to ' + str(filepath.resolve()))
        for x in data:
            print(x, file=f)
        f.close()

def print_cmd(filename, cmd):
    filepath = dirpath / 'results/cmds' / filename
    with filepath.open('w', encoding='utf-8') as f:
        print(cmd, file=f)
        f.close()        


# generate ret[0..n-1] where n users independently select an item from {0..K-1}
# according to Zipf's law with parameter s, and ret[i] is user i's item
def zipfian(n, K, s):
    # calculate the CDF of Zipfian distribution using sentinel
    # so cdf[0] = 0, cdf[K] = 1, and cdf[i] = sum_{j<i} prob[i] more generally
    total = 0
    prob = [0]*K
    for i in range(K):
        prob[i] = (i+1)**(-s)
        total += prob[i]
    prob = [p / total for p in prob]
    cdf = [0]*(K+1)
    for i in range(1, K+1):
        cdf[i] = cdf[i-1] + prob[i-1]

    ret = []
    for _ in range(n):
        while True: # being extra careful in case random() returns 0.0 (unlikely)
            r = random.random()
            if r != 0.0:
                break
        # binary search to find i such that cdf[i-1] < r <= cdf[i]
        i = bisect.bisect_left(cdf, r)
        ret.append(i-1)
    return ret

def zipfian_experiments():
    L = []

    zipf_params = [0.1, 0.5, 1.0, 1.5, 2.0, 3.0]
    for i in range(len(zipf_params)):
        # errors for Zipfian distribution with parameter s
        n = 1000
        K = 22000
        print_data('zipf'+str(zipf_params[i])+'.in', zipfian(n, K, zipf_params[i]))
        zipf = Params(input_file=str((dirpath / 'data' / ('zipf'+str(zipf_params[i])+'.in')).resolve()),
                      output_dir=str((dirpath / 'results').resolve()), num_users=n,
                      run_timing_experiments=False, run_mse_experiments=False, universe_size=K, debug=DEBUG,
                      var_epsilon_min=2.0, var_epsilon_max=4.0, num_var_epsilons=8, include_rappor=False,
                      include_rr=False, include_hpg=False, include_pirappor=False)
        L.append(zipf)
        zipf = Params(input_file=str((dirpath / 'data' / ('zipf'+str(zipf_params[i])+'.in')).resolve()),
                      output_dir=str((dirpath / 'results').resolve()), num_users=n,
                      run_timing_experiments=False, run_variance_experiments=False, universe_size=K, debug=DEBUG,
                      epsilon=5.0, num_trials=300, include_rappor=False, include_rr=False, include_pirappor=False)
        L.append(zipf)
    
    return L

def junta_experiments():
    # user items are equally spread amongst 0,..,junta_size-1
    L = []

    junta_size = [4, 16, 64, 256, 1024, 4096, 8192]
    for i in range(len(junta_size)):
        # errors for Zipfian distribution with parameter s
        n = 1000
        K = 22000
        print_data('junta'+str(junta_size[i])+'.in', [j % junta_size[i] for j in range(n)])
        junta = Params(input_file=str((dirpath / 'data' / ('junta'+str(junta_size[i])+'.in')).resolve()),
                       output_dir=str((dirpath / 'results').resolve()), num_users=n,
                       run_timing_experiments=False, run_mse_experiments=False, universe_size=K, debug=DEBUG,
                       var_epsilon_min=2.0, var_epsilon_max=4.0, num_var_epsilons=8,
                       include_rappor=False,  include_rr=False, include_hpg=False,
                       include_pirappor=False)
        L.append(junta)
        junta = Params(input_file=str((dirpath / 'data' / ('junta'+str(junta_size[i])+'.in')).resolve()),
                       output_dir=str((dirpath / 'results').resolve()), num_users=n,
                       run_timing_experiments=False, run_variance_experiments=False, universe_size=K, debug=DEBUG,
                       epsilon=5.0, num_trials=300, include_rappor=False, include_rr=False, include_pirappor=False)
        L.append(junta)
    
    return L
    
def spike_experiments():
    L = []

    # show that RAPPOR, PI-RAPPOR, PG, and SS have similar error (with SS being best)
    n = 1000
    print_data('spike1.in', [0]*n)
    spike1 = Params(input_file=str((dirpath / 'data' / 'spike1.in').resolve()),
                    output_dir=str((dirpath / 'results').resolve()), num_users=n,
                    run_timing_experiments=False, run_mse_experiments=False,
                    universe_size=1024, debug=DEBUG, var_epsilon_min=2.0, var_epsilon_max=5.5,
                    num_var_epsilons=8, include_hpg=False, include_hr=False, include_rhr=False,
                    include_rr=False)
    L.append(spike1)
    spike1 = Params(input_file=str((dirpath / 'data' / 'spike1.in').resolve()),
                    output_dir=str((dirpath / 'results').resolve()), num_users=n,
                    run_timing_experiments=False, run_variance_experiments=False,
                    epsilon=5.0, universe_size=1024, debug=DEBUG, num_trials=300,
                    include_hpg=False, include_hr=False, include_rhr=False, include_rr=False)
    L.append(spike1)

    # show how much slower RAPPOR is than the others
    n = 1000
    print_data('spike2.in', [0]*n)
    spike2 = Params(input_file=str((dirpath / 'data' / 'spike2.in').resolve()),
                    output_dir=str((dirpath / 'results').resolve()), num_users=n,
                    run_variance_experiments=False, run_mse_experiments=False,
                    debug=DEBUG, universe_min=1024, universe_max=16384,
                    num_universes=20)
#    L.append(spike2)

    # show how much slower PI-RAPPOR is than the others
    n = 1000
    print_data('spike3.in', [0]*n)
    spike3 = Params(input_file=str((dirpath / 'data' / 'spike3.in').resolve()),
                    output_dir=str((dirpath / 'results').resolve()), num_users=n,
                    run_variance_experiments=False, run_mse_experiments=False,
                    debug=DEBUG, universe_min=1024, universe_max=16384,
                    include_rappor=False, num_universes=20)
 #   L.append(spike3)

    # show timing experiments for algorithms other than RAPPOR and PI-RAPPOR
    n = 10000
    print_data('spike4.in', [0]*n)
    spike4 = Params(input_file=str((dirpath / 'data' / 'spike4.in').resolve()),
                    output_dir=str((dirpath / 'results').resolve()), num_users=n,
                    run_variance_experiments=False, run_mse_experiments=False,
                    debug=DEBUG, universe_min=3307948, universe_max=3307948, epsilon=5.0,
                    include_rappor=False, include_ss=False, num_universes=1, hpgq=3)
#    L.append(spike4)

    # show that RR has much worse error than the others
    n = 50000
    print_data('spike5.in', [0]*n)
    spike5 = Params(input_file=str((dirpath / 'data' / 'spike5.in').resolve()),
                    output_dir=str((dirpath / 'results').resolve()), num_users=n,
                    run_timing_experiments=False, run_variance_experiments=False,
                    universe_size=22000, debug=DEBUG, epsilon=5.0, num_trials=300,                
                    include_rappor=False, include_pirappor=False, include_ss=False, hpgq=5)
    L.append(spike5)

    # below is a deprecated experiment with new implementation of HR
    n = 50000
    print_data('spike6.in', [0]*n)
    spike6 = Params(input_file=str((dirpath / 'data' / 'spike6.in').resolve()),
                    output_dir=str((dirpath / 'results').resolve()), num_users=n,
                    run_timing_experiments=False, run_variance_experiments=False,
                    universe_size=22000, debug=DEBUG, epsilon=5.0, num_trials=300,                
                    include_rappor=False, include_pirappor=False, include_rr=False,
                    include_ss=False, hpgq=5)
    L.append(spike6)

    # show that HR and RHR have much worse error than the others, ignoring RR
    n = 10000
    print_data('spike7.in', [0]*n)
    spike7 = Params(input_file=str((dirpath / 'data' / 'spike7.in').resolve()),
                    output_dir=str((dirpath / 'results').resolve()), num_users=n,
                    run_timing_experiments=False, run_variance_experiments=False,
                    universe_size=22000, debug=DEBUG, epsilon=5.0, num_trials=300,                
                    include_rappor=False, include_pirappor=False, include_rr=False,
                    hpgq=5)
    L.append(spike7)

    return L

def run_experiments():
    global dirpath, DEBUG
    dirpath = Path(FLAGS.output_dir)
    DEBUG = FLAGS.debug
    if dirpath.exists() and dirpath.is_dir():
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)
    os.mkdir(dirpath / 'results')
    os.mkdir(dirpath / 'results/cmds')
    os.mkdir(dirpath / 'data')
    L = [] # list of experiments to run, as shell commands
    L += spike_experiments()
    L += zipfian_experiments()
    L += junta_experiments()
    for params in L:
        print("running command: " + str(params))
        base = os.path.basename(params.input_file).rpartition('.')[0]
        ## below assumes each command is only running one type of experiment
        if params.run_mse_experiments: print_cmd(base + '.mse.cmd', str(params))
        elif params.run_variance_experiments: print_cmd(base + '.var.cmd', str(params))
        else: print_cmd(base + '.time.cmd', str(params))
        os.system(str(params))

if __name__ == "__main__":
    app.run(lambda argv: run_experiments())
