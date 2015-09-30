#!/usr/bin/env python

from itertools import tee
import logging

from eden.converter.rna.rnafold import rnafold_to_eden
from eden.modifier.seq import seq_to_seq
from eden.modifier.seq import shuffle_modifier
from eden.util import fit
from eden.util import estimate
from eden.graph import Vectorizer
from eden.converter.fasta import fasta_to_sequence
from eden.util import random_bipartition_iter


logger = logging.getLogger(__name__)


def rfam_url(family_id):
	return 'http://rfam.xfam.org/family/%s/alignment?acc=%s&format=fastau&download=0' % (family_id, family_id)


def binary_classification_dataset_setup(iterable_seq = None, negative_shuffle_ratio = None, shuffle_order = None):

	iter1, iter2 = tee(iterable_seq)
	iterable_graph = rnafold_to_eden(iter1)
	iter3 = seq_to_seq(iter2, modifier = shuffle_modifier, times = negative_shuffle_ratio, order = shuffle_order)
	iterable_graph_neg = rnafold_to_eden(iter3)
	return iterable_graph, iterable_graph_neg


def split_to_train_and_test(rfam_id = None , train_to_test_split_ratio = None):

	iterable = fasta_to_sequence(rfam_url(rfam_id))
	train, test = random_bipartition_iter(iterable, relative_size = train_to_test_split_ratio)
	return train, test


def generate_negatives_and_evaluate(iterable = None , estimator = None, negative_shuffle_ratio = None, shuffle_order = None):

	vectorizer = Vectorizer(complexity=vectorizer_complexity)
	iterable, iterable_neg = binary_classification_dataset_setup(iterable_seq = iterable, negative_shuffle_ratio = negative_shuffle_ratio, shuffle_order = shuffle_order)
	roc, apr =  estimate(iterable, iterable_neg, estimator, vectorizer, n_jobs=-1)
	return roc, apr

	
def generate_negatives_and_fit(iterable = None, negative_shuffle_ratio = None, shuffle_order = None, vectorizer_complexity = None):

	vectorizer = Vectorizer(complexity=vectorizer_complexity)
	iterable, iterable_neg = binary_classification_dataset_setup(iterable_seq = iterable, negative_shuffle_ratio = negative_shuffle_ratio, shuffle_order = shuffle_order)
	model = fit(iterable , iterable_neg , vectorizer , n_jobs=-1 , cv=3 , n_iter_search=1)
	return model


def fit_evaluate(iterable_train = None , iterable_test = None, negative_shuffle_ratio = None, shuffle_order = None, vectorizer_complexity = None):

	estimator = generate_negatives_and_fit(iterable = iterable_train, negative_shuffle_ratio = negative_shuffle_ratio, shuffle_order = shuffle_order, vectorizer_complexity = vectorizer_complexity)
	roc, apr = generate_negatives_and_evaluate(iterable = iterable_test , estimator = estimator, negative_shuffle_ratio = negative_shuffle_ratio, shuffle_order = shuffle_order)
	return roc, apr
	
