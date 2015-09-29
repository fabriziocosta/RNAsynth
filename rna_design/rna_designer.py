#!/usr/bin/env python

import logging
import argparse

from lib.antaRNA import antaRNA_v109


logger = logging.getLogger(__name__)


def get_args():
	"""
	Function for manipulating command-line args.
	Returns a dictionary.
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument('--tGC', '-tgc', type=float, default=0.5, help='')
	parser.add_argument('--colonies', '-c', type=int, default=1, help='')
	parser.add_argument('--name', '-n', type=str, default='antaRNA_', help='')
	parser.add_argument('--alpha', '-a', type=float, default=1.0, help='')
	parser.add_argument('--beta', '-b', type=float, default=1.0, help='')
	parser.add_argument('--evaporation_rate', '-e', type=float, default=0.2, help='')
	parser.add_argument('--struct_correction_term', '-sct', type=float, default=0.5, help='')
	parser.add_argument('--GC_correction_term', '-gct', type=float, default=5.0, help='')
	parser.add_argument('--seq_correction_term', '-qct', type=float, default=1.0, help='')
	parser.add_argument('--degreeOfSequenceInducement', '-dsi', type=int, default=1, help='')
	parser.add_argument('--file_id', '-f', type=str, default='STDOUT', help='')
	parser.add_argument('--verbose', '-v', type=bool, default=False, help='')
	parser.add_argument('--output_verbose', '-o', type=bool, default=False, help='')
	parser.add_argument('--tGCmax', '-tgm', type=float, default=-1.0, help='')
	parser.add_argument('--tGCvar', '-tgv', type=float, default=-1.0, help='')
	parser.add_argument('--termination_convergence', '-t', type=int, default=50, help='')
	parser.add_argument('--convergence_count', '-co', type=int, default=130, help='')
	parser.add_argument('--reset_limit', '-r', type=int, default=5, help='')
	parser.add_argument('--improve', '-i', type=str, default='s', help='')
	parser.add_argument('--temperature', '-tm', type=float, default=37.0, help='')
	parser.add_argument('--paramFile', '-p', type=str, default='', help='')
	parser.add_argument('--return_mod', '-rm', type=bool, default=True, help='')
	parser.add_argument('--seed', '-s', type=str, default='none', help='')
	
	args = parser.parse_args()
	return args


class RNADesign(object):
		
	def __init__(self):
		"""
		DOCUMENTATION
		"""
		self._antaParams = vars(get_args())
		logger.info('Instantiated an RNADesigner object.')
		
	def design(self, dot_bracket_constraint_string, sequence_constraint_string):
		"""
		DOCUMENTATION
		"""
		result = ', '.join(antaRNA_v109.findSequence(dot_bracket_constraint_string, sequence_constraint_string, **self._antaParams))
		sequence = result.split("\n")[2]
		return sequence
	
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
	logger.info('Call to rna_desiner module.')
	designer = RNADesign()
	seq = designer.design("...(((...)))...((...))...","NNNNNNNNNNNNNNNNNNNNNNNNU")
	print seq
