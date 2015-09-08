#!/usr/bin/env python

import logging
from eden import graph
from eden.converter.rna.rnafold import rnafold_to_eden
from eden.converter.fasta import fasta_to_sequence
from itertools import tee, izip
import random
from antaRNA.antaParams import antaParams
from antaRNA import antaRNA_v109
from RNADesign import ConstraintFinder as cf

from sklearn.linear_model import SGDClassifier
from eden.graph import Vectorizer

random_gen_bound = 100000


#logging.basicConfig(level=logging.INFO)


logger = logging.getLogger(__name__)

def pre_process(data, **opts):
    from eden.converter.fasta import fasta_to_sequence
    seqs = fasta_to_sequence(data)
    
    from eden.converter.rna.rnashapes import rnashapes_to_eden
    graphs = rnashapes_to_eden(seqs, shape_type=5, energy_range=35, max_num=opts['max_num'], split_components=True)
	
	return graphs 
	

def generate_antaRNA_sequence(dot_bracket_constraint_string = None, 
							  sequence_constraint_string = None,
							  gc_content = None, 
							  original_header = None, 
							  antaParams = None):
	result = ', '.join(antaRNA_v109.findSequence(dot_bracket_constraint_string,sequence_constraint_string,gc_content,**antaParams.to_dict()))
	header = original_header + '_' + str(random.randrange(random_gen_bound))
	sequence = result.split("\n")[2]
	return header,sequence


class RNASynth(object):
	# self.estimator
	# self.filter_estimator
	# self.vectorizer
	# self.pre_processor
	# self.params
	pass
	
	def init(params,  estimator = SGDClassifier(), vectorizer=Vectorizer(), pre_processor=pre_processor ):
		self.estimator = estimator
		self.vectorizer = vectorizer
		self.pre_processor = pre_processor
		assert param_file is not None, 'ERROR: empty param_file'
		self.params = antaParams(param_file)
		
	#documented
	
	#fit(seqs) train estimators
	#documented
	def fit(self, 
			iterable, 
			iterable_true = None
			negative_shuffle_ratio = None,
			shuffle_order = None):
		"""
		DOCUMENTATION
		"""
		iter1, iter2 = tee(iterable)
		iterable = rnafold_to_eden(iter1)
		iter3 = seq_to_seq(iter2, modifier = shuffle_modifier, times = negative_shuffle_ratio, order = shuffle_order)
		iterable_neg = rnafold_to_eden(iter3)
		
		if iterable_true:
			iter1, iter2 = tee(iterable_true)
			iterable_true = rnafold_to_eden(iter1)
			iter3 = seq_to_seq(iter2, modifier = shuffle_modifier, times = negative_shuffle_ratio, order = shuffle_order)
			iterable_neg_true = rnafold_to_eden(iter3)
			
			iterable = chain(iterable, iterable_neg)
			iterable_neg = chain(iterable_neg , iterable_neg_true)
		
		model = fit(iterable , iterable_neg , self.vectorizer , n_jobs=-1 , cv=3 , n_iter_search=1)
		return model
	
	#_seqs_to_.. private method
	
	def design(self,
			   iterable = None,
			   importance_threshold_sequence_constraint = 0,
			   min_size_connected_component_sequence_constraint = 1,
			   importance_threshold_structure_constraint = 0,
			   min_size_connected_component_structure_constraint = 1,
			   min_size_connected_component_unpaired_structure_constraint = 1,
			   n_synthesized_sequences_per_seed_sequence = 1):
		
		iterable = self.vectorizer.annotate(iterable, estimator = self.estimator)
		iterable = cf.generate_antaRNA_constraints(iterable,importance_threshold_sequence_constraint, 
												   min_size_connected_component_sequence_constraint,
												   importance_threshold_structure_constraint, 
												   min_size_connected_component_structure_constraint, 
												   min_size_connected_component_unpaired_structure_constraint)
		for (dot_bracket,seq_constraint,gc_content,fasta_id)  in iterable:
			for count in range(n_synthesized_sequences_per_seed_sequence):
				result = generate_antaRNA_sequence(dot_bracket, seq_constraint, gc_content, fasta_id , self.params)
				yield result


	def filter(self, iterable = None, threshold = None):
		iterable_sequence, iterable_sequence_for_graphs = tee(iterable)
		graphs = rnafold_to_eden(iterable_sequence_for_graphs)
		predictions = self.vectorizer.predict(graphs , self.filter_estimator)
		for prediction,seq in izip( predictions , iterable_sequence):
			if prediction > threshold:
				yield seq


	#sample(seqs) return synthesized seqs
	#documented
	def sample(self, 
			  iterable = None,
			  importance_threshold_sequence_constraint = None,
			  min_size_connected_component_sequence_constraint = None,
			  importance_threshold_structure_constraint = None,
			  min_size_connected_component_structure_constraint = None,
			  min_size_connected_component_unpaired_structure_constraint = None,
			  n_synthesized_sequences_per_seed_sequence = None,
			  instance_score_threshold = None):
				  
		iterable = self.design(iterable = iterable,
							   importance_threshold_sequence_constraint = importance_threshold_sequence_constraint,
							   min_size_connected_component_sequence_constraint = min_size_connected_component_sequence_constraint,
							   importance_threshold_structure_constraint = importance_threshold_structure_constraint,
							   min_size_connected_component_structure_constraint = min_size_connected_component_structure_constraint,
							   min_size_connected_component_unpaired_structure_constraint = min_size_connected_component_unpaired_structure_constraint,
							   n_synthesized_sequences_per_seed_sequence = n_synthesized_sequences_per_seed_sequence)
		
		iterable = self.filter(iterable = iterable, threshold = instance_score_threshold)
		
		return iterable
			 
###########

def design_RNA(param_file = None , iterable = None, vectorizer = None, estimator = None, \
			importance_threshold_sequence_constraint=0, min_size_connected_component_sequence_constraint=1, importance_threshold_structure_constraint=0, min_size_connected_component_structure_constraint=1, \
			min_size_connected_component_unpaired_structure_constraint=1, n_synthesized_sequences_per_seed_sequence=1):
	"""
	Function for synthesizing RNA sequences.
	Takes as input an iterator over networkx graphs and outputs an iterator 
	over fasta-seq,original-fasta-id pairs.
	Uses antaRNA for sequence synthesis and EDeN for annotating networkx graphs.
	Yields as output [header,sequence].
	"""
	assert param_file is not None, 'ERROR: empty param_file'
	op = antaParams(param_file)

	iterable = vectorizer.annotate(iterable, estimator=estimator)
	iterable = cf.generate_antaRNA_constraints(iterable, importance_threshold_sequence_constraint, min_size_connected_component_sequence_constraint, \
											importance_threshold_structure_constraint, min_size_connected_component_structure_constraint, min_size_connected_component_unpaired_structure_constraint)
	for (dot_bracket,seq_constraint,gc_content,fasta_id)  in iterable:
		for count in range(n_synthesized_sequences_per_seed_sequence):
			result = generate_antaRNA_sequence(dot_bracket, seq_constraint, gc_content, fasta_id ,op)
			yield result


def filter_sequences(iterable = None, vectorizer = None, estimator = None , threshold = 0):
	"""
	Filter. Returns a subset of the iterable with marginal prediction above instance_score_threshold.
	Takes as input a list of fasta sequences. Outputs an iterator over a list of eden sequences.
	"""
	iterable_sequence, iterable_sequence_for_graphs = tee( iterable )
	graphs = rnafold_to_eden( iterable_sequence_for_graphs )
	predictions = vectorizer.predict( graphs , estimator )
	for prediction,seq in izip( predictions , iterable_sequence):
		if prediction > threshold:
			yield seq


def design_filtered_RNA(param_file = None , iterable = None, vectorizer = None, design_estimator = None, filter_estimator = None, \
			importance_threshold_sequence_constraint = 0, min_size_connected_component_sequence_constraint = 1, importance_threshold_structure_constraint = 0, min_size_connected_component_structure_constraint = 1, \
			min_size_connected_component_unpaired_structure_constraint = 1, n_synthesized_sequences_per_seed_sequence = 1 , instance_score_threshold = 0):

	"""
	Synthesize RNA sequences given ..

	Args:
		param_file: 
		iterable: 
		d: the maximal distance size.
		n: the maximal number of clusters used to discretize real label vectors.
		min_r: the minimal radius size.
		min_d: the minimal distance size.
		min_n: the minimal number of clusters used to discretize real label vectors (default 2).
		label_size: the number of discretization steps used in the conversion from real valued labels
		to discrete labels.
		nbits: the number of bits that defines the feature space size:
		|feature space|=2^nbits (default 20).
		normalization: flag to set the resulting feature vector to have unit euclidean norm (default True)
		inner_normalization: flag to set the feature vector for a specific combination of the radius and
		distance size to have unit euclidean norm (default True). When used together with the
		'normalization' flag it will be applied first and then the resulting feature vector
		will be normalized.
		triangular_decomposition: flag to add to each graph the disjoint set of triangles. This
		allows also dense graphs to be processed. Note that runtimes can significantly increase.


	....
	
	Returns:
		descibe

	Takes as input an iterator over networkx graphs and outputs an iterator 
	over fasta-seq,original-fasta-id pairs.
	Uses antaRNA for sequence synthesis and EDeN for annotating networkx graphs.
	Returns as output a fasta list.
	"""
	iterable = design_RNA(param_file = param_file , iterable = iterable , vectorizer = vectorizer, estimator = design_estimator, \
			importance_threshold_sequence_constraint = importance_threshold_sequence_constraint, min_size_connected_component_sequence_constraint = min_size_connected_component_sequence_constraint, \
			importance_threshold_structure_constraint = importance_threshold_structure_constraint, min_size_connected_component_structure_constraint = min_size_connected_component_structure_constraint, \
			min_size_connected_component_unpaired_structure_constraint = min_size_connected_component_unpaired_structure_constraint , n_synthesized_sequences_per_seed_sequence = n_synthesized_sequences_per_seed_sequence)

	iterable = filter_sequences(iterable = iterable, vectorizer = vectorizer, estimator = filter_estimator , threshold = instance_score_threshold)
	return iterable


if __name__ == "__main__":


	logger.info('Call to RNA Design Tools package.')
	
	opts={'antaRNA_param_file':'/home/kohvaeip/RNAsynth/lib/antaRNA/antaRNA.ini' , \
	'importance_threshold_sequence_constraint':-1.1 , 'min_size_connected_component_sequence_constraint':1 , 'importance_threshold_structure_constraint':-0.85 , \
	'min_size_connected_component_structure_constraint':3 , 'min_size_connected_component_unpaired_structure_constraint':3 , 'n_synthesized_sequences_per_seed_sequence':3, 'instance_score_threshold':0}
	
	rfam_id = 'RF00005'
	rfam_url = 'http://rfam.xfam.org/family/%s/alignment?acc=%s&format=fastau&download=0'%(family_id,family_id) 
	iterable = fasta_to_sequence(rfam_url)
	iterable = rnafold_to_eden(iterable)

	#sequences = design_filtered_RNA()
