#!/usr/bin/env python

import argparse
import logging
import time
import sys
from itertools import chain
from itertools import tee
import numpy as np

from eden.util import random_bipartition_iter

from rnasynth.rna_synthesizer import RNASynthesizerInitializer


from util.dataset import split_to_train_and_test
from util.dataset import fit_evaluate


logger = logging.getLogger(__name__)


def performance_evaluation(params, synthesizer=None, iter_train=None, iter_test=None):
    """
    """
    shuffle_order = params['shuffle_order']
    negative_shuffle_ratio = params['negative_shuffle_ratio']
    vectorizer_complexity = params['vectorizer_complexity']

    # Copy training sample iterable for sequence synthesis and producing mixed sample set.
    iter_train, iter_train_syn, iter_seq_true = tee(iter_train, 3)

    # Copy test sample iterable used for evaluation.
    iter_test, iter_test_ = tee(iter_test)

    # Train TrueSamplesModel classifier. Evaluate.
    logger.info('Fit estimator on original data and evaluate the estimator.')
    roc_t, apr_t = fit_evaluate(iterable_train=iter_train, iterable_test=iter_test,
                                negative_shuffle_ratio=negative_shuffle_ratio,
                                shuffle_order=shuffle_order, vectorizer_complexity=vectorizer_complexity)

    # Produce synthesied sequences generator.
    iterable_seq_syn = synthesizer.fit_sample(iter_train_syn)

    # Mix synthesized and true samples.
    iterable_seq_mixed = chain(iterable_seq_syn, iter_seq_true)

    # Train MixedSamplesModel classifier. Evaluate.
    logger.info(
        'Fit estimator on original + sampled data and evaluate the estimator.')
    roc_s, apr_s = fit_evaluate(iterable_train=iterable_seq_mixed, iterable_test=iter_test_,
                                negative_shuffle_ratio=negative_shuffle_ratio,
                                shuffle_order=shuffle_order, vectorizer_complexity=vectorizer_complexity)

    return roc_t, apr_t, roc_s, apr_s


def batch_performance_evaluation(params, synthesizer=None, iter_train=None, iter_test=None, relative_size=None):
    """
    """
    n_experiment_repetitions = params['n_experiment_repetitions']

    start_time = time.time()

    e_roc_t = []
    e_apr_t = []
    e_roc_s = []
    e_apr_s = []

    for epoch in range(n_experiment_repetitions):
        logger.info('-' * 80)
        logger.info('run %d/%d' % (epoch + 1, n_experiment_repetitions))

        # Copy train and test iterables for one run.
        iter_train, iter_train_ = tee(iter_train)
        iter_test, iter_test_ = tee(iter_test)

        # Portion of train and test iterables used in one run.
        iter_train_, x = random_bipartition_iter(
            iter_train_, relative_size=relative_size)
        iter_test_, x = random_bipartition_iter(
            iter_test_, relative_size=relative_size)

        roc_t, apr_t, roc_s, apr_s = performance_evaluation(
            params, synthesizer=synthesizer, iter_train=iter_train_, iter_test=iter_test_)

        # Update experiment performance measures.
        e_roc_t.append(roc_t)
        e_apr_t.append(apr_t)
        e_roc_s.append(roc_s)
        e_apr_s.append(apr_s)

    elapsed_time = time.time() - start_time

    return e_roc_t, e_apr_t, e_roc_s, e_apr_s, elapsed_time


def learning_curve(params):
    """
    """
    rfam_id = params['rfam_id']
    train_to_test_split_ratio = params['train_to_test_split_ratio']

    data_fraction_lower_bound = params['data_fraction_lower_bound']
    data_fraction_upper_bound = params['data_fraction_upper_bound']
    data_fraction_chunks = params['data_fraction_chunks']

    if not(check_data_fractions_integrity(data_fraction_lower_bound, data_fraction_upper_bound,
                                          data_fraction_chunks)):
        sys.stdout.write('Inconsistent data fraction!')
        sys.exit

    data_fractions = list(np.linspace(
        data_fraction_lower_bound, data_fraction_upper_bound, data_fraction_chunks))
    print data_fractions

    logger.info('Starting RNA Synthesis experiment for %s ...' % rfam_id)

    iter_train, iter_test = split_to_train_and_test(
        rfam_id=rfam_id, train_to_test_split_ratio=train_to_test_split_ratio)

    roc_t = []
    roc_s = []
    apr_t = []
    apr_s = []

    for i, data_fraction in enumerate(data_fractions):
        logger.info('=' * 80)
        logger.info('Training on data chunk %d/%d (data fraction: %.1f)' %
                    (i, len(data_fractions), data_fraction))
        iter_train, iter_train_ = tee(iter_train)
        iter_test, iter_test_ = tee(iter_test)

        # Create an instance of the synthesizer object.
        synthesizer = RNASynthesizerInitializer().synthesizer
        mroct, maprt, mrocs, maprs, elapsed_time = batch_performance_evaluation(params,
                                                                                synthesizer=synthesizer,
                                                                                iter_train=iter_train_,
                                                                                iter_test=iter_test_,
                                                                                relative_size=data_fraction)

        logger.info('Performance measures for data fraction: %.1f' %
                    data_fraction)
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


def check_data_fractions_integrity(lower_bound, upper_bound, chuncks):
    if not(0.0 < lower_bound <= 1.0) or not(0.0 <= upper_bound <= 1.0) or (upper_bound <= lower_bound):
        return False
    return True


def get_args():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--rfam_id', '-i', type=str, default='RF00005', help='rfam family ID')
    parser.add_argument('--antaRNA_params', '-p', type=str,
                        default='./antaRNA.ini', help='antaRNA initialization file')
    parser.add_argument('--n_experiment_repetitions', '-j',
                        type=int, default=10, help='runs per experiment')
    parser.add_argument('--train_to_test_split_ratio', '-r',
                        type=float, default=0.2, help='train-to-test size ratio')
    parser.add_argument(
        '--shuffle_order', '-s', type=int, default=2, help='shuffle order')
    parser.add_argument('--negative_shuffle_ratio', '-ns', type=int,
                        default=2, help='number of negative samples generated per sample')
    parser.add_argument(
        '--vectorizer_complexity', '-v', type=int, default=2, help='eden Vectorizer complexity')
    parser.add_argument('--data_fraction_lower_bound', '-dfl',
                        type=float, default=0.1,
                        help='lower bound for generating the list of sample data fractions')
    parser.add_argument(
        '--data_fraction_upper_bound', '-dfu', type=float, default=1.0,
        help='upper bound for generating the list of sample data fractions')
    parser.add_argument(
        '--data_fraction_chunks', '-dfc', type=int, default=10,
        help='number of chunks in the list of sample data fractions')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    """
    """
    logging.basicConfig(level=logging.INFO)
    logger.info('Call to performance_evaluation module.')

    params = vars(get_args())

    roc_t, roc_s, apr_t, apr_s, data_fractions = learning_curve(params)
