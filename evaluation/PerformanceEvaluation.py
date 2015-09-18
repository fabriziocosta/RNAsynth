#!/usr/bin/env python

import argparse
import logging
import time
import sys
from itertools import chain
from itertools import tee

from eden.util import fit
from eden.util import estimate
from eden.graph import Vectorizer
from eden.converter.fasta import fasta_to_sequence
from eden.converter.rna.rnafold import rnafold_to_eden
from eden.modifier.seq import seq_to_seq
from eden.modifier.seq import shuffle_modifier
from eden.util import random_bipartition_iter

from rna_design.rna_synthesizer import RNASynth

from util.dataset import rfam_url
from util.dataset import binary_classification_dataset_setup
from util.dataset import split_to_train_and_test
from util.dataset import generate_negatives_and_evaluate
from util.dataset import generate_negatives_and_fit


logger = logging.getLogger(__name__)


def learning_curve(params, iter_train = None, iter_test = None, relative_size=None):
    """
    DOCUMENTATION
    
    Executes one epoch of n runs.
    Each run does the experiment on "relative_size" of training samples.
    Outputs four lists containing the overall ROC and APR performance measures for true and mixed training samples.
    Also outputs the overall elapsed time on n runs.
    """
	
    log_file = params['log_file']
    antaRNA_param_file = params['antaRNA_params']
    n_experiment_repetitions = params['n_experiment_repetitions']
    shuffle_order = params['shuffle_order']
    negative_shuffle_ratio = params['negative_shuffle_ratio']
    vectorizer_complexity = params['vectorizer_complexity']

    start_time = time.time()

    e_roc_t = []
    e_apr_t = []
    e_roc_s = []
    e_apr_s = []


    for epoch in range(n_experiment_repetitions):
        logger.info('-' * 80)
        logger.info('run %d/%d' % (epoch+1, n_experiment_repetitions))
        
        # copy train and test iterables for one run
        iter_train, iter_train_ = tee(iter_train)
        iter_test, iter_test_ = tee(iter_test)

        # portion of train and test iterables used in one run
        iter_train_ , x = random_bipartition_iter(iter_train_, relative_size=relative_size)
        iter_test_ , x = random_bipartition_iter(iter_test_, relative_size=relative_size)
        
        # train iterable used for sequence synthesis and producing mixed samples set.
        iter_train_, iter_train_syn, iter_true = tee(iter_train_, 3)
        
        # test iterable used for evaluation.
        iter_test_, iter_test__ = tee(iter_test_)
        
        # Train TrueSamplesModel classifier.
        logger.info('Fit estimator on original data')
        estimator_true = generate_negatives_and_fit(iterable = iter_train_, negative_shuffle_ratio = negative_shuffle_ratio,
											   shuffle_order = shuffle_order, vectorizer_complexity = vectorizer_complexity)
        
        # evaluate the true samples model -> output ROC and APR performance measures
        logger.info('Evaluate true samples estimator:')
        roc_t, apr_t = generate_negatives_and_evaluate(iterable = iter_test_ , estimator = estimator_true,
													   negative_shuffle_ratio = negative_shuffle_ratio, shuffle_order = shuffle_order)

        # Create an RNASynth object.
        synthesizer = RNASynth(params)
        
        # Produce synthesied sequences generator.
        iterable_seq_syn = synthesizer.sample(iterable_seq = iter_train_syn)
        
        # Mix synthesized and true samples.
        iterable_seq_mixed = chain(iterable_seq_syn, iter_seq_true)

        logger.info('Fit estimator on original + sampled data')
        # Train MixedSamplesModel classifier.
        estimator_mixed = generate_negatives_and_fit(iterable = iterable_seq_mixed, negative_shuffle_ratio = negative_shuffle_ratio,
													 shuffle_order = shuffle_order, vectorizer_complexity = vectorizer_complexity)

        # evaluate the mixed samples model -> output ROC and APR performance measures
        logger.info('Evaluate mixed samples estimator:')
        roc_s, apr_s = generate_negatives_and_evaluate(iterable = iter_test__ , estimator = estimator_mixed,
													   negative_shuffle_ratio = negative_shuffle_ratio, shuffle_order = shuffle_order)
        
        # Update experiment performance measures:
        e_roc_t.append(roc_t)
        e_apr_t.append(apr_t)
        e_roc_s.append(roc_s)
        e_apr_s.append(apr_s)

    elapsed_time = time.time() - start_time

    return e_roc_t, e_apr_t, e_roc_s, e_apr_s, elapsed_time


def get_args():
    """
    Function for manipulating command-line args.
    Returns a dictionary.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--rfam_id', '-i', type=str, default='RF00005', help='rfam family ID')
    parser.add_argument('--log_file', '-l', type=str, default='~/Synthesis.log', help='experiment log file')
    parser.add_argument('--antaRNA_params', '-p', type=str, default='./antaRNA.ini', help='antaRNA initialization file')
    parser.add_argument('--importance_threshold_sequence_constraint', '-a', type=int, default=0, help='nucleotide selection threshold')
    parser.add_argument('--min_size_connected_component_sequence_constraint', '-b', type=int, default=1, help='nucleotide minimum adjacency')
    parser.add_argument('--importance_threshold_structure_constraint', '-c', type=int, default=0, help='basepair selection threshold')
    parser.add_argument('--min_size_connected_component_structure_constraint', '-d', type=int, default=1, help='basepairs minimum adjacency')
    parser.add_argument('--min_size_connected_component_unpaired_structure_constraint', '-e', type=int, default=1, help='unpaired nucleotides minimum adjacency')
    parser.add_argument('--n_synthesized_sequences_per_seed_sequence', '-n', type=int, default=1, help='number of synthesized sequences per constraint')
    parser.add_argument('--instance_score_threshold', '-f', type=int, default=0, help='filtering threshold')
    parser.add_argument('--n_experiment_repetitions', '-j', type=int, default=10, help='runs per experiment')
    parser.add_argument('--train_to_test_split_ratio', '-r', type=float, default=0.2, help='train-to-test size ratio')
    parser.add_argument('--shuffle_order', '-s', type=int, default=2, help='shuffle order')
    parser.add_argument('--negative_shuffle_ratio', '-ns', type=int, default=2, help='number of negative samples generated per sample')
    parser.add_argument('--vectorizer_complexity', '-v', type=int, default=2, help='????????')
    #parser.add_argument('--data_fractions', '-df', type=list, default=[50,100], help='????????')
    args = parser.parse_args()
    return args


def compute_learning_curves(params):
    """
    Main body of RNASynthesis experiment.
    """
    
    # Experiment parameters
    rfam_id = params['rfam_id']
    log_file = params['log_file']
    train_to_test_split_ratio = params['train_to_test_split_ratio']
    data_fractions = params['data_fractions']
    
    logger.info('Starting RNA Synthesis experiment for %s ...' % rfam_id)
    
    iter_train, iter_test = split_to_train_and_test(rfam_id = rfam_id , train_to_test_split_ratio = train_to_test_split_ratio)
	
    roc_t = []
    roc_s = []
    apr_t = []
    apr_s = []


    for i,data_fraction in enumerate(data_fractions):
        logger.info('=' * 80)
        logger.info('Training on data chunk %d/%d (data fraction: %.1f)' % (i,len(data_fractions),data_fraction))
        iter_train, iter_train_ = tee(iter_train)
        iter_test, iter_test_ = tee(iter_test)

        mroct, maprt, mrocs, maprs, elapsed_time = learning_curve(params, iter_train = iter_train_, iter_test = iter_test_, 
																  relative_size = data_fraction)
  
        logger.info('Performance measures for data fraction: %.1f' % data_fraction)
        logger.info('ROC for True samples:')
        logger.info('%s' % mroct)
        logger.info('ROC for Mixed samples:')
        logger.info('%s' % mrocs)
        logger.info('APR for True samples:')
        logger.info('%s' % maprt)
        logger.info('APR for Mixed samples:')
        logger.info('%s' % maprs)
        logger.info('Elapsed time: %s' % elapsed_time)
        roc_t.append(mroct)
        roc_s.append(mrocs)
        apr_t.append(maprt)
        apr_s.append(maprs)
    return roc_t, roc_s, apr_t, apr_s, data_fractions


if __name__ == "__main__":
    """
    Main body of the experiment.
    Runs on a specific rfam family.
    Starts from a portion of the training samples and increases the portion per iteration.
    Runs an epoch on each portion.
    Saves the results in individual files labeld with the portion.
    """
    logging.basicConfig(level=logging.INFO)
    logger.info('Call to experiment module.')

    params = vars(get_args())

    roc_t, roc_s, apr_t, apr_s = compute_learning_curves(params)
