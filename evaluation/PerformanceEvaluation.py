#!/usr/bin/env python

from itertools import chain
from itertools import tee
import time
import sys
from eden.util import fit
from eden.util import estimate
from eden.graph import Vectorizer
from eden.converter.fasta import fasta_to_sequence
from eden.converter.rna.rnafold import rnafold_to_eden
from eden.util import random_bipartition_iter
from eden.modifier.seq import seq_to_seq
from eden.modifier.seq import shuffle_modifier
from eden.util import random_bipartition_iter
from RNADesign import RNADesignTools as rdt

import argparse
import logging
logger = logging.getLogger(__name__)


def rfam_url(family_id):
    return 'http://rfam.xfam.org/family/%s/alignment?acc=%s&format=fastau&download=0' % (family_id, family_id)


def produce_batch_RNA(graphs=None, vectorizer=None, estimator=None,
                      antaRNA_param_file=None, importance_threshold_sequence_constraint=0, min_size_connected_component_sequence_constraint=1,
                      importance_threshold_structure_constraint=0, min_size_connected_component_structure_constraint=1, min_size_connected_component_unpaired_structure_constraint=1,
                      n_synthesized_sequences_per_seed_sequence=1, instance_score_threshold=0):

    graphs, graphs_counter = tee(graphs)

    # Design the batch of new sequences.
    fasta_iterable = rdt.design_filtered_RNA(param_file=antaRNA_param_file,
                                             iterable=graphs,
                                             vectorizer=vectorizer,
                                             design_estimator=estimator,
                                             filter_estimator=estimator,
                                             importance_threshold_sequence_constraint=importance_threshold_sequence_constraint,
                                             min_size_connected_component_sequence_constraint=min_size_connected_component_sequence_constraint,
                                             importance_threshold_structure_constraint=importance_threshold_structure_constraint,
                                             min_size_connected_component_structure_constraint=min_size_connected_component_structure_constraint,
                                             min_size_connected_component_unpaired_structure_constraint=min_size_connected_component_unpaired_structure_constraint,
                                             n_synthesized_sequences_per_seed_sequence=n_synthesized_sequences_per_seed_sequence,
                                             instance_score_threshold=instance_score_threshold)

    for i, seq in enumerate(fasta_iterable):
        yield seq


def run_epoch(n_experiment_repetitions=1, relative_size=None, antaRNA_param_file=None, performance_log_file=None,
              graphs_pos_test=None, graphs_neg_test=None, graphs_pos_train=None, graphs_neg_train=None,
              importance_threshold_sequence_constraint=0, min_size_connected_component_sequence_constraint=1,
              importance_threshold_structure_constraint=0, min_size_connected_component_structure_constraint=1, min_size_connected_component_unpaired_structure_constraint=1,
              n_synthesized_sequences_per_seed_sequence=1, instance_score_threshold=0, vectorizer_complexity=2,
              negative_shuffle_ratio=2):
    """
    Executes one epoch of n runs.
    Each run does the experiment on "relative_size" of training samples.
    Outputs four lists containing the overall ROC and APR performance measures for true and mixed training samples.
    Also outputs the overall elapsed time on n runs.
    """

    opts = {'importance_threshold_sequence_constraint': importance_threshold_sequence_constraint, 
			'min_size_connected_component_sequence_constraint': min_size_connected_component_sequence_constraint,
            'importance_threshold_structure_constraint': importance_threshold_structure_constraint, 'min_size_connected_component_structure_constraint': min_size_connected_component_structure_constraint,
            'min_size_connected_component_unpaired_structure_constraint': min_size_connected_component_unpaired_structure_constraint, 'n_synthesized_sequences_per_seed_sequence': n_synthesized_sequences_per_seed_sequence,
            'instance_score_threshold': instance_score_threshold, 'antaRNA_param_file': antaRNA_param_file}

    start_time = time.time()

    measure_list_ROCT = []
    measure_list_APRT = []
    measure_list_ROCS = []
    measure_list_APRS = []

    # Instantiate the vectorizer.
    vectorizer = Vectorizer(complexity=vectorizer_complexity)

    for epoch in range(n_experiment_repetitions):
        logger.info('-' * 80)
        logger.info('run %d/%d' % (epoch+1, n_experiment_repetitions))
        graphs_pos_test, graphs_pos_test_epoch_1, graphs_pos_test_epoch_2 = tee(graphs_pos_test, 3)
        graphs_neg_test, graphs_neg_test_epoch_1, graphs_neg_test_epoch_2 = tee(graphs_neg_test, 3)

        # Get a fresh copy of test and train graphs for the new epoch.
        graphs_pos_train, graphs_pos_train_ = tee(graphs_pos_train)
        graphs_neg_train, graphs_neg_train_ = tee(graphs_neg_train)

        # Take n% of test and train graph sets. n is the epoch counter.
        graphs_pos_train_epoch, iterable_not_used = random_bipartition_iter(graphs_pos_train_, relative_size=relative_size)

        graphs_neg_train_epoch, iterable_not_used = random_bipartition_iter(graphs_neg_train_, relative_size=relative_size)

        # Make copies of each for the whole epoch run.
        graphs_pos_train_epoch_1, graphs_pos_train_epoch_2, graphs_pos_train_epoch_3 = tee(graphs_pos_train_epoch, 3)

        graphs_neg_train_epoch_1, graphs_neg_train_epoch_2 = tee(graphs_neg_train_epoch)

        # Train TrueSamplesModel classifier.
        logger.info('Fit estimator on original data')
        estimator = fit(graphs_pos_train_epoch_1, graphs_neg_train_epoch_1, vectorizer, n_jobs=-1, cv=3, n_iter_search=1)

        # Test the test set against TrueSamplesModel -> output performance
        logger.info('Evaluate estimator:')
        roc_t, apr_t = estimate(graphs_pos_test_epoch_1, graphs_neg_test_epoch_1, estimator, vectorizer, n_jobs=-1)

        # Design the batch of new sequences.
        designed_iterable = produce_batch_RNA(graphs=graphs_pos_train_epoch_2, vectorizer=vectorizer, estimator=estimator, **opts)

        # Make copies of designed_iterable.
        designed_iterable, designed_iterable2 = tee(designed_iterable)

        # Build the graphs iterator over the designed sequences.
        graphs_synthesized = rnafold_to_eden(designed_iterable)

        # Produce graph iterator over negative samples for synthesized data.
        iterable_neg = seq_to_seq(designed_iterable2, modifier=shuffle_modifier, times=negative_shuffle_ratio, order=2)
        graphs_neg_synthesized = rnafold_to_eden(iterable_neg)

        # Mix the sample with "true sequences".
        graphs_mixed_pos = chain(graphs_pos_train_epoch_3, graphs_synthesized)
        graphs_mixed_neg = chain(graphs_neg_train_epoch_2, graphs_neg_synthesized)

        logger.info('Fit estimator on original + sampled data')
        # Train MixedSamplesModel classifier.
        estimator2 = fit(graphs_mixed_pos, graphs_mixed_neg, vectorizer, n_jobs=-1, cv=3, n_iter_search=1)

        # Test the test set against SynthesizedSamplesModel -> output performance
        logger.info('Evaluate estimator:')
        roc_s, apr_s = estimate(graphs_pos_test_epoch_2, graphs_neg_test_epoch_2, estimator2, vectorizer, n_jobs=-1)

        # Additional sample counts:
        measure_list_ROCT.append(roc_t)
        measure_list_APRT.append(apr_t)
        measure_list_ROCS.append(roc_s)
        measure_list_APRS.append(apr_s)

    elapsed_time = time.time() - start_time

    return measure_list_ROCT, measure_list_APRT, measure_list_ROCS, measure_list_APRS, elapsed_time


def get_args():
    """
    Function for manipulating command-line args.
    Returns a dictionary.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--rfam_id', '-i', type=str, default='RF00005', help='rfam family ID')
    # --log_file is not currently used.
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
    args = parser.parse_args()
    return args


def compute_learning_curves(params):
    """
    Main body of RNASynthesis experiment.
    """
    logger.info('Starting RNA Synthesis experiment for %s ...' % params['rfam_id'])

    iterable_pos = fasta_to_sequence(rfam_url(params['rfam_id']))
    iterable_pos, iterable_pos_ = tee(iterable_pos)
    iterable_neg = seq_to_seq(iterable_pos_, modifier=shuffle_modifier, times=params['negative_shuffle_ratio'], order=2)

    # Positive sample graphs.
    graphs_pos = rnafold_to_eden(iterable_pos)
    # Negative sample graphs.
    graphs_neg = rnafold_to_eden(iterable_neg)
    # Split samples to train and test iterables.
    graphs_pos_train, graphs_pos_test = random_bipartition_iter(graphs_pos, relative_size=params['train_to_test_split_ratio'])
    graphs_neg_train, graphs_neg_test = random_bipartition_iter(graphs_neg, relative_size=params['train_to_test_split_ratio'])

    exp_roc_t = []
    exp_roc_s = []
    exp_apr_t = []
    exp_apr_s = []

    for i,data_fraction in enumerate(params["data_fractions"]):
        logger.info('=' * 80)
        logger.info('Training on data chunk %d/%d (data fraction: %.3f)' % (i,len(params["data_fractions"]),data_fraction))
        graphs_pos_train, graphs_pos_train_c = tee(graphs_pos_train)
        graphs_pos_test, graphs_pos_test_c = tee(graphs_pos_test)
        graphs_neg_train, graphs_neg_train_c = tee(graphs_neg_train)
        graphs_neg_test, graphs_neg_test_c = tee(graphs_neg_test)

        mroct, maprt, mrocs, maprs, elapsed_time = run_epoch(graphs_pos_test=graphs_pos_test_c,
                                                             graphs_neg_test=graphs_neg_test_c,
                                                             graphs_pos_train=graphs_pos_train_c,
                                                             graphs_neg_train=graphs_neg_train_c,
                                                             importance_threshold_sequence_constraint=params['importance_threshold_sequence_constraint'],
                                                             min_size_connected_component_sequence_constraint=params['min_size_connected_component_sequence_constraint'],
                                                             importance_threshold_structure_constraint=params['importance_threshold_structure_constraint'],
                                                             min_size_connected_component_structure_constraint=params['min_size_connected_component_structure_constraint'],
                                                             min_size_connected_component_unpaired_structure_constraint=params['min_size_connected_component_unpaired_structure_constraint'],
                                                             n_synthesized_sequences_per_seed_sequence=params['n_synthesized_sequences_per_seed_sequence'],
                                                             instance_score_threshold=params['instance_score_threshold'],
                                                             n_experiment_repetitions=params['n_experiment_repetitions'],
                                                             relative_size=data_fraction,
                                                             antaRNA_param_file=params['antaRNA_params'],
                                                             vectorizer_complexity=params['vectorizer_complexity'],
                                                             negative_shuffle_ratio=params['negative_shuffle_ratio'])

        logger.info('Performance measures for data fraction: %.3f' % data_fraction)
        logger.info('ROC for True samples:')
        logger.info('%s' % mroct)
        logger.info('ROC for Mixed samples:')
        logger.info('%s' % mrocs)
        logger.info('APR for True samples:')
        logger.info('%s' % maprt)
        logger.info('APR for Mixed samples:')
        logger.info('%s' % maprs)
        logger.info('Elapsed time: %s' % elapsed_time)
        exp_roc_t.append(mroct)
        exp_roc_s.append(mrocs)
        exp_apr_t.append(maprt)
        exp_apr_s.append(maprs)
    return exp_roc_t, exp_roc_s, exp_apr_t, exp_apr_s, params["data_fractions"]


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

    roc_t, roc_s, apr_t, apr_s = experiment(params)
