#!/usr/bin/env python

#
#
#          genSeqErrorModel.py
#          Computes sequencing error model for gen_reads.py
#
#         
#          Usage: python genSeqErrorModel.py -i input_reads.fq -o path/to/output_name.p
#
#
# Python 3 ready

# Testing Git

import numpy as np
import pandas as pd
import argparse
import sys
import pickle
import matplotlib.pyplot as mpl
import pathlib
import pysam
from functools import reduce

# enables import from neighboring package

# sys.path.append(str(pathlib.Path(__file__).resolve().parents[1]))
# from source.probability import DiscreteDistribution


sys.path.append('/Users/ralhazmy/Documents/NEAT/NEAT/source')
#sys.path.insert(1, '/Users/ralhazmy/Documents/NEAT/NEAT')
from probability import DiscreteDistribution

#from NEAT.source.probability import DiscreteDistribution

from bisect import bisect_left

def take_closest(bins, quality):
    """
    Assumes bins is sorted. Returns closest value to quality.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(bins, quality)
    if pos == 0:
        return bins[0]
    if pos == len(bins):
        return bins[-1]
    before = bins[pos - 1]
    after = bins[pos]
    if after - quality < quality - before:
        return after
    else:
        return before


def parse_file(input_file, quality_scores, real_q, off_q, max_reads, n_samp):
    # Takes a gzip or sam file and returns the simulation's average error rate,
    print('reading ' + input_file + '...')
    is_aligned = False
    lines_to_read = 0
    try:
        if input_file[-4:] == '.bam' or input_file[-4:] == '.sam':
            print('detected aligned file....')
            stats = pysam.idxstats(input_file).strip().split('\n')
            lines_to_read = reduce(lambda x, y: x + y, [eval('+'.join(l.rstrip('\n').split('\t')[2:])) for l in stats])
            f = pysam.AlignmentFile(input_file)
            is_aligned = True
        else:
            print('detected fastq file....')
            with pysam.FastxFile(input_file) as f:
                for _ in f:
                    lines_to_read += 1
            f = pysam.FastxFile(input_file)
    except FileNotFoundError:
        print("Check input file. Must be fastq, gzipped fastq, or bam/sam file.")
        sys.exit(1)

    actual_readlen = 0
    #q_dict = {}
    current_line = 0
    quarters = lines_to_read // 4

    #A list of n dictionaries, n being the length of a sequence // Put outside of this loop
    error_model = {
        'quality_scores': quality_scores, #a list of q scores
        'quality_score_probabilities': [], # number of bins * length of sequence data frame to hold the quality score probabilities
        'quality_offset': off_q,
        'avg_error': {},
        'error_parameters': [
                            # sequencing substitution transition probabilities
                            [[0., 0.4918, 0.3377, 0.1705], [0.5238, 0., 0.2661, 0.2101], [0.3754, 0.2355, 0., 0.3890], [0.2505, 0.2552, 0.4942, 0.]],
                            # if a sequencing error occurs, what are the odds it's an indel?
                            0.01,
                            # sequencing indel error length distribution
                            [0.999, 0.001],
                            [1, 2],
                            # if a sequencing indel error occurs, what are the odds it's an insertion as opposed to a deletion?
                            0.4,
                            # if a sequencing insertion error occurs, what's the probability of it being an A, C, G, T...
                            [0.25, 0.25, 0.25, 0.25]
                            ]
    }


    if is_aligned:
        g = f.fetch()
    else:
        g = f

    obtained_read_length = False #Used to get read length, the read length will = the fist read length.
    temp_q_count = 0

    for read in g:
        if is_aligned:
            qualities_to_check = read.query_alignment_qualities
        else:
            qualities_to_check = read.get_quality_array()
        if actual_readlen == 0: #read length = the length of the first read
            actual_readlen = len(qualities_to_check) - 1
            print('assuming read length is uniform...')
            print('detected read length (from first read found):', actual_readlen)
            
        # check if read length is more than 0 and if we have the read length already**    
        if actual_readlen > 0 and not obtained_read_length: 
            temp_q_count = np.zeros((actual_readlen, len(quality_scores)))
            error_model['quality_score_probabilities'] = np.zeros((actual_readlen, len(quality_scores)), dtype= float)
            obtained_read_length = True

        # sanity-check readlengths
        if len(qualities_to_check) - 1 != actual_readlen:
            print('skipping read with unexpected length...')
            continue
        
        for i in range(0, actual_readlen):
            q = qualities_to_check[i] #The qualities of each base
            bin = take_closest(quality_scores, q)
            bin_index = quality_scores.index(bin)
            temp_q_count[i][bin_index] +=1

        #loading
        current_line += 1
        if current_line % quarters == 0:
            print(f'{(current_line/lines_to_read)*100:.0f}%')
        if 0 < max_reads <= current_line:
            break
    f.close()

    #porbability calculator
    base_index = 0
    for base in temp_q_count:#for every dictionary(scores of a single base in the reads)
        total = sum(base)
        bin_index = 0
        for score in base:
            error_model['quality_score_probabilities'][base_index][bin_index] = score / total
            bin_index += 1
        base_index += 1

    # A Discrete distribution of the quality_score probabilities
    Discretes = []
    for base in error_model['quality_score_probabilities']:
        Discretes.append(DiscreteDistribution(error_model['quality_scores'], base))

    # A counter for the Discrete distribution results run for n_samp times
    count_dict = {}
    for q in error_model['quality_scores']:
        count_dict[q] = 0
    #example: {0: 0, 12: 0, 24: 0, 36: 0}
    lines_to_sample = len(range(1, n_samp + 1))
    samp_quarters = lines_to_sample // 4 #divide the reads into 1/4s for the loading bar
    for samp in range(1, n_samp + 1):
        if samp % samp_quarters == 0:
            print(f'{(samp/lines_to_sample)*100:.0f}%') #loading bar
        for i in range(actual_readlen):
            my_q = Discretes[i].sample()
            count_dict[my_q] += 1

    print(count_dict)

    #Calculates the average error rate
    tot_bases = float(sum(count_dict.values()))
    avg_err = 0.
    for k in sorted(count_dict.keys()):
        eVal = 10. ** (-k / 10.)
        #print (k, eVal, count_dict[k])
        avg_err += eVal * (count_dict[k] / tot_bases)
    print('AVG ERROR RATE:', avg_err)

    return error_model


def main():
    parser = argparse.ArgumentParser(description='genSeqErrorModel.py')
    parser.add_argument('-i', type=str, required=True, metavar='<str>', help="* input_read1.fq (.gz) / input_read1.sam")
    parser.add_argument('-o', type=str, required=True, metavar='<str>', help="* output.p")
    # parser.add_argument('-i2', type=str, required=False, metavar='<str>', default=None,
    #                     help="input_read2.fq (.gz)")
    # parser.add_argument('-p', type=str, required=False, metavar='<str>', default=None, help="input_alignment.pileup")
    parser.add_argument('-q', type=int, required=False, metavar='<int>', default=33, help="quality score offset [33]")
    parser.add_argument('-Q', type=int, required=False, metavar='<int>', default=[0, 12, 24, 36], help="List of quality score bins [0, 12, 24, 36]")
    parser.add_argument('-n', type=int, required=False, metavar='<int>', default=-1,
                        help="maximum number of reads to process [all]")
    parser.add_argument('-s', type=int, required=False, metavar='<int>', default=1000000,
                        help="number of simulation iterations [1000000]")
    parser.add_argument('--plot', required=False, action='store_true', default=False,
                        help='perform some optional plotting')
    args = parser.parse_args()

    (infile, outfile, off_q, max_q, max_reads, n_samp) = (args.i, args.o, args.q, args.Q, args.n, args.s)
    (infile2, pile_up) = (args.i2, args.p)

    real_q = max_q + 1

    plot_stuff = args.plot

    q_scores = range(real_q)
    if infile2 is None:
        (init_q, prob_q, avg_err) = parse_file(infile, real_q, off_q, max_reads, n_samp, plot_stuff)
    else:
        (init_q, prob_q, avg_err1) = parse_file(infile, real_q, off_q, max_reads, n_samp, plot_stuff)
        (init_q2, prob_q2, avg_err2) = parse_file(infile2, real_q, off_q, max_reads, n_samp, plot_stuff)
        avg_err = (avg_err1 + avg_err2) / 2.


    #
    #	embed some default sequencing error parameters if no pileup is provided
    #
    if pile_up == None:
        
        print('Using default sequencing error parameters...')

        # sequencing substitution transition probabilities
        sse_prob = [[0., 0.4918, 0.3377, 0.1705],
                    [0.5238, 0., 0.2661, 0.2101],
                    [0.3754, 0.2355, 0., 0.3890],
                    [0.2505, 0.2552, 0.4942, 0.]]
        # if a sequencing error occurs, what are the odds it's an indel?
        sie_rate = 0.01
        # sequencing indel error length distribution
        sie_prob = [0.999, 0.001]
        sie_val = [1, 2]
        # if a sequencing indel error occurs, what are the odds it's an insertion as opposed to a deletion?
        sie_ins_freq = 0.4
        # if a sequencing insertion error occurs, what's the probability of it being an A, C, G, T...
        sie_ins_nucl = [0.25, 0.25, 0.25, 0.25]

    #
    #	otherwise we need to parse a pileup and compute statistics!
    #
    else:
        print('\nPileup parsing coming soon!\n')
        exit(1)

    err_params = [sse_prob, sie_rate, sie_prob, sie_val, sie_ins_freq, sie_ins_nucl]





    #
    #	finally, let's save our output model
    #
    outfile = pathlib.Path(outfile).with_suffix(".p")
    print('saving model...')
    if infile2 is None:
        #pickle.dump({quality_scores:[], quality_scores_probabilities:[]})
        pickle.dump([init_q, prob_q, q_scores, off_q, avg_err, err_params], open(outfile, 'wb'))
    else:
        pickle.dump([init_q, prob_q, init_q2, prob_q2, q_scores, off_q, avg_err, err_params], open(outfile, 'wb'))

#pickle

if __name__ == '__main__':
    main()
