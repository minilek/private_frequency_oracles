#!/usr/bin/env python3

import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from sys import stdin

import argparse
import os
from pathlib import Path
import shutil
import sys
from absl import app,flags

from pfoparams import Params

FLAGS = flags.FLAGS

flags.DEFINE_string("results_dir", "pfotemp", "directory where results live")

colors = ["b", "g", "r", "c", "mediumslateblue", "y", "dimgrey", "magenta", "mediumspringgreen", "gainsboro", "rosybrown", "chocolate", "aqua", "mediumslateblue", "saddlebrown", "darkorchid", "khaki"]
alg_to_color = {}
color_idx = 0
    

def extract_base_filename(s):
    s = os.path.basename(s)
    return s[:s.rfind('.')]

def get_params_from_file(filepath, ext):
    global cmds_path
    cmdpath = cmds_path / (extract_base_filename(str(filepath.resolve())) + '.' + ext + '.cmd')
    with cmdpath.open('r', encoding='utf-8') as f:
        return Params.string_to_param(f.readline().rstrip())

def gen_var_fig(filepath):
    global colors, alg_to_color, color_idx
    print("opening " + str(filepath.resolve()), file=sys.stderr)
    algos = []
    variances = []
    my_colors = []
    algo_name_to_idx = {}

    params = get_params_from_file(filepath, 'var')
    test_type = extract_base_filename(str(filepath.resolve()))
    if test_type.startswith('spike'):
        test_type = 'spike'

    max_variance = -1

    with filepath.open('r', encoding='utf-8') as f:
        for line in f:
            a = line.split(",")
            alg, eps, var = a[0], float(a[1]), float(a[2])
            max_variance = max(max_variance, var)
            if alg not in algo_name_to_idx:
                algo_name_to_idx[alg] = len(algos)
                if alg not in alg_to_color:
                    alg_to_color[alg] = colors[color_idx]
                    color_idx += 1
                algos.append(alg)
                my_colors.append(alg_to_color[alg])
                variances.append({})
            idx = algo_name_to_idx[alg]
            variances[idx][eps] = var

    for a in algos:
        if set(variances[algo_name_to_idx[a]].keys()) != set(variances[0].keys()):
            print('Will not create figure for ' + filepath.resolve() + ' since not all algos were ran with same values of epsilon.')
            return

    all_epsilons = list(variances[0].keys())
    all_epsilons.sort()

    fig=plt.figure()
    plt.xlabel('epsilon')
    plt.ylabel('MSE')

    fig.text(0.5, .95, test_type+',k='+str(params.universe_size)+',n='+str(params.num_users), horizontalalignment='center', verticalalignment='top')

    plots = []
    for i in range(len(algos)):
        plots.append([])
        for j in range(len(all_epsilons)):
            plots[i].append(variances[i][all_epsilons[j]])
        plt.plot(all_epsilons, plots[i], '-', label=algos[i], color=my_colors[i])

    lgnd = plt.legend(loc='upper right')

    for color,text in zip(my_colors, lgnd.get_texts()):
        text.set_color(color)
    
    plt.savefig(str((figs_path / (extract_base_filename(str(filepath.resolve())) + '.var.png')).resolve()), facecolor='floralwhite', edgecolor='none')

def gen_mse_or_max_fig(filepath, mse_type=True):
    global colors, alg_to_color, color_idx
    print("opening " + str(filepath.resolve()), file=sys.stderr)
    algos = []
    errors = []
    my_colors = []
    algo_name_to_idx = {}

    params = get_params_from_file(filepath, 'mse')
    test_type = extract_base_filename(str(filepath.resolve()))
    if test_type.startswith('spike'):
        test_type = 'spike'

    max_error = -1

    with filepath.open('r', encoding='utf-8') as f:
        for line in f:
            a = line.split(",")
            alg, err = a[0], float(a[1])
            max_error = max(max_error, err)
            if alg not in algo_name_to_idx:
                algo_name_to_idx[alg] = len(algos)
                if alg not in alg_to_color:
                    alg_to_color[alg] = colors[color_idx]
                    color_idx += 1
                algos.append(alg)
                my_colors.append(alg_to_color[alg])
                errors.append([])
            idx = algo_name_to_idx[alg]
            errors[idx].append(err)

    for i in range(len(errors)):
        errors[i].sort()

    # plot CDFs to within .01 granularity on x-axis
    N = len(errors[0])
    dx = 1.0/N
    X  = np.arange(0, 100.0, dx*100)

    fig=plt.figure()
    plt.xlabel('percentile')
    if mse_type:
        plt.ylabel('MSE')
    else:
        plt.ylabel('max error')

    fig.text(0.5, .95, test_type+',k='+str(params.universe_size)+',n='+str(params.num_users)+',eps='+str(params.epsilon), horizontalalignment='center', verticalalignment='top')

    plots = []
    for i in range(len(algos)):
        plots.append([])
        for j in range(N):
            plots[i].append(errors[i][int(1.0*j/N * len(errors[i]))])
        plt.plot(X, plots[i], 'o', label=algos[i], markersize=2, color=my_colors[i])

    lgnd = plt.legend(loc='upper left')

    for color,text in zip(my_colors, lgnd.get_texts()):
        text.set_color(color)

    suffix = '.mse.png'
    if not mse_type:
        suffix = '.max.png'
    plt.savefig(str((figs_path / (extract_base_filename(str(filepath.resolve())) + suffix)).resolve()), facecolor='floralwhite', edgecolor='none')

def gen_time_fig(filepath):
    global colors, alg_to_color, color_idx
    print("opening " + str(filepath.resolve()), file=sys.stderr)
    algos = []
    runtimes = []
    my_colors = []
    algo_name_to_idx = {}

    params = get_params_from_file(filepath, 'time')
    test_type = extract_base_filename(str(filepath.resolve()))
    if test_type.startswith('spike'):
        test_type = 'spike'

    max_time = -1

    with filepath.open('r', encoding='utf-8') as f:
        for line in f:
            a = line.split(",")
            alg, k, time = a[0], int(a[1]), int(float(a[2])) # float() first in case # is so big it uses scientific notation
            max_time = max(max_time, time)
            if alg not in algo_name_to_idx:
                algo_name_to_idx[alg] = len(algos)
                if alg not in alg_to_color:
                    alg_to_color[alg] = colors[color_idx]
                    color_idx += 1
                algos.append(alg)
                my_colors.append(alg_to_color[alg])
                runtimes.append({})
            idx = algo_name_to_idx[alg]
            runtimes[idx][k] = time

    for a in algos:
        if set(runtimes[algo_name_to_idx[a]].keys()) != set(runtimes[0].keys()):
            print('Will not create figure for ' + filepath.resolve() + ' since not all algos were ran with same values of k.')
            return

    all_universe_sizes = list(runtimes[0].keys())
    all_universe_sizes.sort()

    fig=plt.figure()
    plt.xlabel('universe size')
    plt.ylabel('runtime')

    fig.text(0.5, .95, test_type+',n='+str(params.num_users)+',eps='+str(params.epsilon), horizontalalignment='center', verticalalignment='top')

    plots = []
    for i in range(len(algos)):
        plots.append([])
        for j in range(len(all_universe_sizes)):
            plots[i].append(runtimes[i][all_universe_sizes[j]])
        plt.plot(all_universe_sizes, plots[i], '-', label=algos[i], color=my_colors[i])

    lgnd = plt.legend(loc='upper left')

    for color,text in zip(my_colors, lgnd.get_texts()):
        text.set_color(color)
    
    plt.savefig(str((figs_path / (extract_base_filename(str(filepath.resolve())) + '.time.png')).resolve()), facecolor='floralwhite', edgecolor='none')

def generate_figures():
    global results_path, cmds_path, figs_path

    SMALL_SIZE = 8
    MEDIUM_SIZE = 18
    BIGGER_SIZE = 24

    plt.rcParams["figure.figsize"] = (11.5,7)
    plt.rc('font', size=BIGGER_SIZE)         # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    
    results_path = Path(FLAGS.results_dir)
    cmds_path = results_path / 'cmds'
    figs_path = results_path.parents[0] / 'figures'
    if figs_path.exists():
        shutil.rmtree(figs_path)
    os.mkdir(figs_path)
    for f in results_path.iterdir():
        if not f.is_dir():
            pathstr = str(f.resolve())
            if pathstr.endswith('.time'):
                gen_time_fig(f)
            elif pathstr.endswith('.mse') or pathstr.endswith('.max'):
                gen_mse_or_max_fig(f, pathstr.endswith('.mse'))
            elif pathstr.endswith('.var'):
                gen_var_fig(f)
            else:
                print('unknown result filetype ' + str(f.resolve()))

if __name__ == "__main__":
    app.run(lambda argv: generate_figures())
