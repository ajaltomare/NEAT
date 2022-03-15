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
sys.path.append(str(pathlib.Path(__file__).resolve().parents[1]))

from source.probability import DiscreteDistribution

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


def parse_file(input_file, real_q, off_q, max_reads, n_samp, plot_stuff):
    init_smooth = 0.
    prob_smooth = 0.

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
    q_dict = {}
    current_line = 0
    quarters = lines_to_read // 4

    input_string = input('Enter the bins wanted, spereated by a space:\n Example:0 12 24 36')
    if False: #change later 
        sys.exit('Error: Invalid input.\n exiting...')

    #Makes a dictionary with the user's bins as keys
    inputs_list= input_string.split()
    quality_bins = {} #list of user inputs      *** 
    for bins in inputs_list:
        quality_bins[bins] = 0

    #A list of n dictionaries, n being the length of a sequence // Put outside of this loop
    error_model = {
        'q_scores': [quality_bins], #a list of q scores dictionary for each base
        'q_score_probabilities': np.zeros(actual_readlen, columns = inputs_list, dtype= float), # number of bins * length of sequence data frame to hold the quality score probabilities
        'quality_offset': {off_q}, #(off_q)??
        'avg_error': {},
        'error_parameters': {}#hard coded
    }


    if is_aligned:
        g = f.fetch()
    else:
        g = f

    obtained_read_length = False #Used to get read length, the read length will = the fist read length.

    for read in g:
        if is_aligned:
            qualities_to_check = read.query_alignment_qualities
        else:
            qualities_to_check = read.get_quality_array()
        if actual_readlen == 0:
            actual_readlen = len(qualities_to_check) - 1
            print('assuming read length is uniform...')
            print('detected read length (from first read found):', actual_readlen)
            
        # check if read length is more than 0 and if we have the read length already**    
        if actual_readlen > 0 and not obtained_read_length: 
            error_model['q_scores'] = [quality_bins] * (actual_readlen)
            obtained_read_length = True

        # sanity-check readlengths
        if len(qualities_to_check) - 1 != actual_readlen:
            print('skipping read with unexpected length...')
            continue
        
    
        for i in range(0, actual_readlen-1):
            q = qualities_to_check[i] #The qualities of each base
            q_dict[q] = True #??
            bin = take_closest(quality_bins, q)
            error_model['q_scores'][i][bin] +=1


        current_line += 1
        if current_line % quarters == 0:
            print(f'{(current_line/lines_to_read)*100:.0f}%')
        if 0 < max_reads <= current_line:
            break

    f.close()


    #porbability calculator
    pdIndex = 0
    for dict in error_model['q_scores']:#for every dictionary(scores of a single base in the reads)
        total = sum(dict.values())
        if not total == actual_readlen: # total should = the read length
            print('error')# among other things
        for key, value in dict.items(): # items() returns a set-like view backed by the dict
            error_model['q_score_probabilities'][pdIndex][key] = value / total
        pdIndex += 1


    # some sanity checking again...
    q_range = [min(q_dict.keys()), max(q_dict.keys())]
    if q_range[0] < 0:
        print('\nError: Read in Q-scores below 0\n')
        exit(1)
    if q_range[1] > real_q:
        print('\nError: Read in Q-scores above specified maximum:', q_range[1], '>', real_q, '\n')
        exit(1)


    print('estimating average error rate via simulation...')
    q_scores = range(real_q)
    # print (len(init_q), len(init_q[0]))
    # print (len(prob_q), len(prob_q[1]), len(prob_q[1][0]))

    init_dist_by_pos = [DiscreteDistribution(init_q[i], q_scores) for i in range(len(init_q))]
    prob_dist_by_pos_by_prev_q = [None]
    for i in range(1, len(init_q)):
        prob_dist_by_pos_by_prev_q.append([])
        for j in range(len(init_q[0])):
            if np.sum(prob_q[i][j]) <= 0.:  # if we don't have sufficient data for a transition, use the previous qscore
                prob_dist_by_pos_by_prev_q[-1].append(DiscreteDistribution([1], [q_scores[j]], degenerate_val=q_scores[j]))
            else:
                prob_dist_by_pos_by_prev_q[-1].append(DiscreteDistribution(prob_q[i][j], q_scores))

#average error
    count_dict = {} 
    for q in q_scores: #q_scores is a range(real_q), which is an argument.
        count_dict[q] = 0
    lines_to_sample = len(range(1, n_samp + 1)) #n_samp is an argument, number of reads??
    samp_quarters = lines_to_sample // 4 #divide the reads into 1/4s ??
    for samp in range(1, n_samp + 1):
        if samp % samp_quarters == 0:
            print(f'{(samp/lines_to_sample)*100:.0f}%') #loading bar or something like that
        my_q = init_dist_by_pos[0].sample() #was a discrete distribution
        count_dict[my_q] += 1
        for i in range(1, len(init_q)):
            my_q = prob_dist_by_pos_by_prev_q[i][my_q].sample() #works with prevoius position, are we still doing that?
            count_dict[my_q] += 1

    tot_bases = float(sum(count_dict.values()))
    avg_err = 0.
    for k in sorted(count_dict.keys()):
        eVal = 10. ** (-k / 10.)
        # print k, eVal, count_dict[k]
        avg_err += eVal * (count_dict[k] / tot_bases)
    print('AVG ERROR RATE:', avg_err)

    quality_score_probability_matrix = error_model['q_score_probabilities']
    error_model['q_score_probabilities'] = quality_score_probability_matrix.apply(
       lambda row: DiscreteDistribution(error_model['quality_scores'], row), axis=1).to_numpy()


    return init_q, prob_q, avg_err


def main():
    parser = argparse.ArgumentParser(description='genSeqErrorModel.py')
    parser.add_argument('-i', type=str, required=True, metavar='<str>', help="* input_read1.fq (.gz) / input_read1.sam")
    parser.add_argument('-o', type=str, required=True, metavar='<str>', help="* output.p")
    # parser.add_argument('-i2', type=str, required=False, metavar='<str>', default=None,
    #                     help="input_read2.fq (.gz)")
    # parser.add_argument('-p', type=str, required=False, metavar='<str>', default=None, help="input_alignment.pileup")
    # parser.add_argument('-q', type=int, required=False, metavar='<int>', default=33, help="quality score offset [33]")
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
