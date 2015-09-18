#!/usr/bin/env python

#sys

#libr

#local

from itertools import tee
from itertools import izip
import logging

from sklearn.linear_model import SGDClassifier

from eden.converter.rna.rnafold import rnafold_to_eden
from eden.converter.fasta import fasta_to_sequence
from eden.converter.rna.rnashapes import rnashapes_to_eden
from eden.util import fit
from eden.modifier.seq import seq_to_seq
from eden.modifier.seq import shuffle_modifier
from eden.graph import Vectorizer

from antaRNA.antaParams import antaParams
from antaRNA import antaRNA_v109
from rna_design import constraint_extractor as ce

from util.dataset import rfam_url
from util.dataset import binary_classification_dataset_setup


logger = logging.getLogger(__name__)


def pre_process(iterable_seq, **opts):
	"""
	DOCUMENTATION
	"""
	graphs = rnashapes_to_eden(seqs, shape_type=5, energy_range=35, max_num=opts['max_num'], split_components=True)
	return graphs 


def generate_antaRNA_sequence(dot_bracket_constraint_string = None, 
							  sequence_constraint_string = None,
							  gc_content = None, 
							  original_header = None,
							  index = None, 
							  antaParams = None):
	result = ', '.join(antaRNA_v109.findSequence(dot_bracket_constraint_string,sequence_constraint_string,gc_content,**antaParams.to_dict()))
	header = original_header + '_' + str(index)
	sequence = result.split("\n")[2]
	return header,sequence


class RNASynth(object):
		
	def __init__(self, params,  estimator = SGDClassifier(), vectorizer = Vectorizer(), pre_processor = pre_process):
		
		self._antaRNA_param_file = params['antaRNA_params']
		self._importance_threshold_sequence_constraint = params['importance_threshold_sequence_constraint']
		self._min_size_connected_component_sequence_constraint = params['min_size_connected_component_sequence_constraint']
		self._importance_threshold_structure_constraint = params['importance_threshold_structure_constraint']
		self._min_size_connected_component_structure_constraint = params['min_size_connected_component_structure_constraint']
		self._min_size_connected_component_unpaired_structure_constraint = params['min_size_connected_component_unpaired_structure_constraint']
		self._n_synthesized_sequences_per_seed_sequence = params['n_synthesized_sequences_per_seed_sequence']
		self._instance_score_threshold = params['instance_score_threshold']
		self._train_to_test_split_ratio = params['train_to_test_split_ratio']
		self._shuffle_order = params['shuffle_order']
		self._negative_shuffle_ratio = params['negative_shuffle_ratio']
		self._vectorizer_complexity = params['vectorizer_complexity']
		self._max_num_graphs_per_seq = params['max_num_graphs_per_seq']
		
		self.estimator = estimator
		self.vectorizer = vectorizer
		self.pre_processor = pre_processor
		
		assert self._antaRNA_param_file is not None, 'ERROR: empty param_file'
		self._antaParams = antaParams(self._antaRNA_param_file)
		logger.info('Instantiated an RNASynthesizer object.')
		
	def __fit(self, iterable_seq):
		"""
		DOCUMENTATION
		"""
		iterable_graph, iterable_graph_neg = binary_classification_dataset_setup(iterable_seq = iterable_seq, negative_shuffle_ratio = self._negative_shuffle_ratio, 
																	 shuffle_order = self._shuffle_order)
		model = fit(iterable_graph, iterable_graph_neg , self.vectorizer , n_jobs=-1 , cv=3 , n_iter_search=1)
		return model
	
	def __design(self, iterable_graph):
		"""
		DOCUMENTATION
		"""
		iterable_graph = self.vectorizer.annotate(iterable_graph, estimator = self.estimator)
		iterable = ce.generate_antaRNA_constraints(iterable_graph, self._importance_threshold_sequence_constraint, 
												   self._min_size_connected_component_sequence_constraint,
												   self._importance_threshold_structure_constraint, 
												   self._min_size_connected_component_structure_constraint, 
												   self._min_size_connected_component_unpaired_structure_constraint)
		for (dot_bracket,seq_constraint,gc_content,fasta_id)  in iterable:
			for count in range(self._n_synthesized_sequences_per_seed_sequence):
				result = generate_antaRNA_sequence(dot_bracket, seq_constraint, gc_content, fasta_id , count, self._antaParams)
				#result = self.designer.design(dot_bracket, seq_constraint, gc_content, fasta_id , count, self._antaParams)
				yield result

	def __filter(self, iterable_seq):
		"""
		DOCUMENTATION
		"""
		iter1, iter2 = tee(iterable_seq)
		iterable_graph = rnafold_to_eden(iter1)
		predictions = self.vectorizer.predict(iterable_graph , self.estimator)
		
		for prediction,seq in izip(predictions, iter2):
			if prediction > self._instance_score_threshold:
				yield seq

	def sample(self, iterable_seq):
		"""
		DOCUMENTATION
		"""	
		iterable_seq, iterable_seq_ = tee(iterable_seq)
		self.estimator = self.__fit(iterable_seq)
		iterable_graph = rnafold_to_eden(iterable_seq_)
		iterable_seq = self.__design(iterable_graph)
		iterable_seq = self.__filter(iterable_seq)
		
		return iterable_seq
			 

	"""
	Synthesize RNA sequences given ..

	Args:
		antaRNA_params: parameter setting file for antaRNA wrappper.
		importance_threshold_sequence_constraint: classification score threshold for identifying important nucleotides in a sequence.
		min_size_connected_component_sequence_constraint: minimum number of adjacent important nucleotides which can form a sequence constraint.
		importance_threshold_structure_constraint: classification score threshold for labeling important basepairs in a secondary structure.
		min_size_connected_component_structure_constraint: minimum number of adjacent basepairs which can form a secondary structure constraint. 
		This constraint is in the form of dot-bracket notation
		min_size_connected_component_unpaired_structure_constraint: 
		n_synthesized_sequences_per_seed_sequence: option for setting the number of synthesized sequences per constraint.
		instance_score_threshold: predicted score threshold for filtering synthesized sequences.
		train_to_test_split_ratio: ratio for splitting the sample dataset into train and test datasets.
		shuffle_order: eden.modifier.seq.seq_to_seq parameter. 
		negative_shuffle_ratio: number of negative sample sequences generated for each positive sample.
		vectorizer_complexity: eden.graph.Vectorizer parameter.
		max_num_graphs_per_seq: eden.converter.rna.rnashapes.rnashapes_to_eden parameter.

	....
	
	Returns:
		descibe

	Takes as input an iterator over networkx graphs and outputs an iterator 
	over fasta-seq,original-fasta-id pairs.
	Uses antaRNA for sequence synthesis and EDeN for annotating networkx graphs.
	Returns as output a fasta list.
	"""

	
if __name__ == "__main__":

	logging.basicConfig(level=logging.INFO)
	logger.info('Call to RNASynthesizer module.')
	
	params = {'antaRNA_params':'./antaRNA.ini', 'importance_threshold_sequence_constraint':0, 
		  'min_size_connected_component_sequence_constraint':1, 'importance_threshold_structure_constraint':0,
		  'min_size_connected_component_structure_constraint':1, 'min_size_connected_component_unpaired_structure_constraint':1,
		  'n_synthesized_sequences_per_seed_sequence':2, 'instance_score_threshold':0,
		  'train_to_test_split_ratio':.2, 'shuffle_order':2, 'negative_shuffle_ratio':2,
		  'vectorizer_complexity':2, 'max_num_graphs_per_seq':3}

	rfam_id = 'RF01685'
	iterable_seq = fasta_to_sequence('http://rfam.xfam.org/family/%s/alignment?acc=%s&format=fastau&download=0' % (rfam_id, rfam_id))
	synthesizer = RNASynth(params)
	iter_seq = synthesizer.sample(iterable_seq)
	for item in iter_seq:
		print item
