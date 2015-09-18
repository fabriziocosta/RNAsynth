#!/usr/bin/env python

from itertools import tee
import logging

from eden.converter.rna.rnafold import rnafold_to_eden
from eden.modifier.seq import seq_to_seq
from eden.modifier.seq import shuffle_modifier


logger = logging.getLogger(__name__)


def rfam_url(family_id):
	return 'http://rfam.xfam.org/family/%s/alignment?acc=%s&format=fastau&download=0' % (family_id, family_id)


def binary_classification_dataset_setup(iterable_seq = None, negative_shuffle_ratio = None, shuffle_order = None):
	"""
	DOCUMENTATION
	"""
	iter1, iter2 = tee(iterable_seq)
	iterable_graph = rnafold_to_eden(iter1)
	iter3 = seq_to_seq(iter2, modifier = shuffle_modifier, times = negative_shuffle_ratio, order = shuffle_order)
	iterable_graph_neg = rnafold_to_eden(iter3)
	return iterable_graph, iterable_graph_neg
