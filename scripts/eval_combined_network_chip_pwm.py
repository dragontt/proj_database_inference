#!/usr/bin/python

"""
Combine mutliple methods evaluated by chip support and pwm support.
The figures are in form of cumulative predictions.

eval_method: cumulative or binned
"""

import sys
import argparse
import glob
import os.path
import numpy
import matplotlib.pyplot as plt

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Combine chip and pwm evaluations")
    parser.add_argument('-e', '-eval_method', dest='eval_method', type=str)
    parser.add_argument('-n', '-fn_np', dest='fn_np', type=str)
    parser.add_argument('-b', '-fn_bart', dest='fn_bart', type=str)
    parser.add_argument('-i', '-dir_infer', dest='dir_infer', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    # list of methods
    methods = ['dbl_flt_pid_cutoff_40', 'dbl_flt_pid_knn_5', 'dbl_flt_pid_knn_pct10', 'quantile_normalization_standard', 'scertf_pwm']

    # figure setup
    offset = 2
    methods.extend(['netprophet', 'holstege_bart_np'])
    # methods.extend(['netprophet', 'fire_np'])
    colors = ['magenta', 'cyan', 'green', 'yellow', 'red', 'black', 'blue']
    x_ticks = ['4k', '8k', '12k', '16k', '20k', '24k', '28k', '32k', '36k', '40k']

    # compute chip and pwm support ratios of each method
    parsed.dir_infer = check_dir(parsed.dir_infer)
    eval_chip = [None] * len(methods)
    eval_pwm = [None] * len(methods)

    [eval_chip[len(methods)-offset], eval_pwm[len(methods)-offset]] = \
        parse_ratio(parsed.fn_np, parsed.eval_method)
    [eval_chip[len(methods)-offset+1], eval_pwm[len(methods)-offset+1]] = \
        parse_ratio(parsed.fn_bart, parsed.eval_method)

    for i in range(len(methods)-offset):
        # print methods[i]
        [eval_chip[i], eval_pwm[i]] = parse_ratio(parsed.dir_infer + methods[i] \
            + "/chip.bp.np.set.sizes.top4to40k.txt", parsed.eval_method)

    # plot figures
    plt.figure(num=None, figsize=(15,8), dpi=80)
    plt.subplot(1,2,1)
    for i in range(len(eval_chip)):
        plt.plot(eval_chip[i], color=colors[i], label=methods[i])
    plt.xticks(range(len(eval_chip[0])), x_ticks)
    plt.xlabel('Predictions grouped by rank')
    plt.ylabel('Interactions supported by ChIP')
    plt.xlim(-1, len(eval_chip[0]))
    plt.ylim(0, 0.35)
    plt.legend(loc="upper right")

    plt.subplot(1,2,2)
    for i in range(len(eval_pwm)):
        plt.plot(eval_pwm[i], color=colors[i], label=methods[i])
    plt.xticks(range(len(eval_pwm[0])), x_ticks)
    plt.xlabel('Predictions grouped by rank')
    plt.ylabel('Interactions supported by PWM')
    plt.xlim(-1, len(eval_pwm[0]))
    plt.ylim(0, 0.35)
    plt.legend(loc="upper right")

    plt.show()

def check_dir(fd):
    if not fd.endswith('/'):
        fd += '/'
    return fd

def parse_ratio(fn, method):
    lines = open(fn, "r").readlines()
    chip = [0] * (len(lines)/3)
    pwm = [0] * (len(lines)/3)
    if method == "cumulative":
        for i in range(len(lines)/3):
            line = lines[i*3].split()
            chip[i] = float(line[5])/float(line[2])
            pwm[i] = float(line[4])/float(line[2])
    elif method == "binned":
        for i in range(len(lines)/3):
            line = lines[i*3].split()
            if i == 0:
                chip[i] = float(line[5])/float(line[2])
                pwm[i] = float(line[4])/float(line[2])
            else:
                chip[i] = (float(line[5]) - float(prevline[5]))/(float(line[2]) - float(prevline[2]))
                pwm[i] = (float(line[4]) - float(prevline[4]))/(float(line[2]) - float(prevline[2]))
            prevline = line    
    return  [chip, pwm]

if __name__ == "__main__":
    main(sys.argv)
