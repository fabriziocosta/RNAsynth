#!/usr/bin/env python

import logging
from antaParams import antaParams
import antaRNA_v109
import sys
from eden import graph
from eden.converter.rna.rnafold import rnafold_to_eden
from eden.converter.fasta import fasta_to_sequence
from eden.modifier.seq import seq_to_seq,  shuffle_modifier
from itertools import tee
import random
import time
import ConstraintFinder as cf

random_gen_bound = 100000


def design_RNA(param_file=None , iterable=None, vectorizer=None, estimator=None, \
			nt_importance_threshold=0, nmin_important_nt_adjaceny=1, bp_importance_threshold=0, nmin_important_bp_adjaceny=1, \
			nmin_unpaired_nt_adjacency=1, multi_sequence_size=1):

	"""
	Function for synthesizing RNA sequences.
	Takes as input an iterator over networkx graphs and outputs an iterator 
	over fasta-seq,original-fasta-id pairs.
	Uses antaRNA for sequence synthesis and EDeN for annotating networkx graphs.
	Returns as output a fasta list.
	"""
	op = antaParams(param_file)

	iterable = vectorizer.annotate(iterable, estimator=estimator)
	iterable = cf.generate_antaRNA_constraints(iterable, nt_importance_threshold, nmin_important_nt_adjaceny, \
											bp_importance_threshold, nmin_important_bp_adjaceny, nmin_unpaired_nt_adjacency)
	output_list = []
	for (dot_bracket,seq_constraint,gc_content,fasta_id)  in iterable:
		for count in range(multi_sequence_size):
			result = generate_antaRNA_sequence(dot_bracket, seq_constraint, gc_content, fasta_id ,op)
			output_list = output_list + result
	return output_list

 
def sequence_to_fasta(sequences):
	fasta_list = []
	for sequence in sequences:
		fasta_list.append(sequence[0])
		fasta_list.append(sequence[1])
	return fasta_list


def generate_antaRNA_sequence(dot_bracket_constraint_string=None,sequence_constraint_string=None,\
							gc_content=None,original_fasta_id=None,antaParams=None):
	result = ', '.join(antaRNA_v109.findSequence(dot_bracket_constraint_string,sequence_constraint_string,gc_content,**antaParams.to_dict()))
	header = original_fasta_id + '_' + str(random.randrange(random_gen_bound))
	sequence = result.split("\n")[2]
	return [header,sequence]


def filter_sequences(iterable=None, vectorizer=None, filtering_estimator=None , filtering_threshold = 0):
	"""
	Filter. Returns a subset of the iterable with marginal prediction above filtering_threshold.
	Takes as input a list of fasta sequences. Outputs an iterator over a list of fasta sequences.
	"""
	iterable_sequence = fasta_to_sequence( iterable )
	iterable_sequence, iterable_sequence_for_graphs, iterable_sequence_for_headers = tee( iterable_sequence, 3 )
	graphs = rnafold_to_eden( iterable_sequence_for_graphs )

	predictions = vectorizer.predict( graphs , filtering_estimator )
	prediction_list = [prediction for prediction in predictions]
	fasta_list = [seq_line for seq_line in iterable_sequence_for_headers]
	candidate_list = zip( prediction_list , fasta_list)
	filtered_sequences = []
	for candidate in candidate_list:
		if candidate[0] > filtering_threshold:
			filtered_sequences.append(candidate[1])
	#todo: message in case of no survivors
	return sequence_to_fasta(filtered_sequences)


def design_batch_RNA(param_file=None , iterable=None, vectorizer=None, design_estimator=None, filter_estimator=None , \
					nt_importance_threshold=0, nmin_important_nt_adjaceny=1, bp_importance_threshold=0, nmin_important_bp_adjaceny=1, \
					nmin_unpaired_nt_adjacency=1, batch_size=1, multi_sequence_size=1, filtering_threshold = 0, sequence_pool=None):
	"""
	Function for synthesizing RNA sequences with filtering.
	Takes as input an iterator over EDeN graphs and outputs 
	an iterator over a fasta list.
	Filters generated sequences using filter_estimator.
	"""
	logging.basicConfig(level=logging.INFO)
	logger = logging.getLogger(__name__)
	
	timestamp_start = time.time()
	op = antaParams(param_file)

	iterable = vectorizer.annotate(iterable, estimator=design_estimator)
	iterable_constraints = cf.generate_antaRNA_constraints(iterable, nt_importance_threshold, nmin_important_nt_adjaceny, \
														bp_importance_threshold, nmin_important_bp_adjaceny, nmin_unpaired_nt_adjacency)
	batch_sequence_list = []
	output_list = []
	duplicate_counter = 0

	while (len(batch_sequence_list) < batch_size*2) and (duplicate_counter < len(sequence_pool)) :
		fasta_list = []
		iterable_constraints, iterable = tee(iterable_constraints)
		#sys.stdout.write('new while iteration with duplicate counter %d\n' %duplicate_counter)
		for (dot_bracket,seq_constraint,gc_content,fasta_id)  in iterable:
			#sys.stdout.write(dot_bracket)
			#sys.stdout.write(seq_constraint)
			#sys.stdout.write(fasta_id)
			for count in range(multi_sequence_size):
				result = generate_antaRNA_sequence(dot_bracket, seq_constraint, gc_content, fasta_id ,op)
				if (not(result[1] in sequence_pool)):
					fasta_list = fasta_list + result
					sequence_pool.append(result[1])
				else:
					sys.stdout.write('duplicate found: %s \n' %result[1])
					duplicate_counter += 1
			if duplicate_counter >= len(sequence_pool):
				logger.info('%d duplicates produced.' %duplicate_counter)
				logger.info('Duplicates exceeded the number of existing samples. Quitting batch synthesis...')
				break
			if 	len(batch_sequence_list) >= batch_size*2:
				break
				
		if len(fasta_list) > 0:
			batch_sequence_list = batch_sequence_list + filter_sequences(iterable=fasta_list, vectorizer=vectorizer, filtering_estimator=filter_estimator , filtering_threshold = 0)

	timestamp_elapsed = time.time() - timestamp_start

	if len(batch_sequence_list) >= batch_size*2 :
		output_list = batch_sequence_list[:batch_size*2]
	else:
		output_list = batch_sequence_list
		logger.info('Design Batch RNA produced %d rather than %d sequences.' %(len(output_list),batch_size*2))
	
	return output_list, timestamp_elapsed


if __name__ == "__main__":
	
	logging.basicConfig(level=logging.INFO)
	logger = logging.getLogger(__name__)
	
	logger.info('Call to RNADesignTools package.')
