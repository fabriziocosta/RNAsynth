#!/usr/bin/env python

from itertools import chain
from itertools import tee
import time
import logging
import sys
from eden.util import fit
from eden.util import estimate
from eden.graph import Vectorizer
import RNADesignTools as rdt
from eden.converter.fasta import fasta_to_sequence
from eden.converter.rna.rnafold import rnafold_to_eden
from eden.util import random_bipartition_iter
from eden.modifier.seq import seq_to_seq
from eden.modifier.seq import shuffle_modifier
from eden.util import random_bipartition_iter
import argparse


def rfam_url(family_id):
	return 'http://rfam.xfam.org/family/%s/alignment?acc=%s&format=fastau&download=0'%(family_id,family_id)


def produce_batch_RNA(graphs = None, vectorizer = None, estimator = None, \
					antaRNA_param_file = None , nt_importance_threshold=0, nmin_important_nt_adjaceny=1, \
					bp_importance_threshold=0, nmin_important_bp_adjaceny=1, nmin_unpaired_nt_adjacency=1, \
					synthesized_batch_proportion=1, multi_sequence_size=1, filtering_threshold = 0):

	graphs , graphs_counter = tee(graphs)

	# Calculate the batch size proportional to the size of the epoch's training samples.
	train_size = 0
	for graph in graphs_counter:
		train_size += 1
	#sys.stdout.write('Number of training samples in this epoch: %d\n' %train_size)
	batch_size = train_size * synthesized_batch_proportion
	# Design the batch of new sequences.
	fasta_iterable = rdt.design_filtered_RNA(param_file = antaRNA_param_file , iterable = graphs, \
											vectorizer = vectorizer, design_estimator = estimator, \
											filter_estimator = estimator , nt_importance_threshold = nt_importance_threshold, \
											nmin_important_nt_adjaceny = nmin_important_nt_adjaceny, bp_importance_threshold = bp_importance_threshold, \
											nmin_important_bp_adjaceny = nmin_important_bp_adjaceny, nmin_unpaired_nt_adjacency = nmin_unpaired_nt_adjacency, \
											multi_sequence_size = multi_sequence_size , filtering_threshold = filtering_threshold)

	for i, seq in enumerate(fasta_iterable):
		if i > batch_size * 2:
			break
		yield seq


def run_epoch(runs = None , relative_size = None, antaRNA_param_file = None , performance_log_file = None , \
		graphs_pos_test = None , graphs_neg_test = None , graphs_pos_train = None , graphs_neg_train = None , \
		nt_importance_threshold=0, nmin_important_nt_adjaceny=1, \
		bp_importance_threshold=0, nmin_important_bp_adjaceny=1, nmin_unpaired_nt_adjacency=1, \
		synthesized_batch_proportion=1, multi_sequence_size=1, filtering_threshold = 0 ):
	"""
	Executes one epoch of n runs.
	Each run does the experiment on "relative_size" of training samples.
	Outputs four lists containing the overall ROC and APR performance measures for true and mixed training samples.
	Also outputs the overall elapsed time on n runs.
	"""

	opts = {'nt_importance_threshold':nt_importance_threshold, 'nmin_important_nt_adjaceny':nmin_important_nt_adjaceny,\
			'bp_importance_threshold':bp_importance_threshold, 'nmin_important_bp_adjaceny':nmin_important_bp_adjaceny, \
			'nmin_unpaired_nt_adjacency':nmin_unpaired_nt_adjacency, 'synthesized_batch_proportion':synthesized_batch_proportion, \
			'multi_sequence_size':multi_sequence_size, 'filtering_threshold':filtering_threshold, 'antaRNA_param_file':antaRNA_param_file }

	start_time = time.time()

	measure_list_ROCT = []
	measure_list_APRT = []
	measure_list_ROCS = []
	measure_list_APRS = []

	for epoch in range(1,runs + 1):

		graphs_pos_test, graphs_pos_test_epoch_1, graphs_pos_test_epoch_2 = tee(graphs_pos_test,3)
		graphs_neg_test, graphs_neg_test_epoch_1, graphs_neg_test_epoch_2 = tee(graphs_neg_test,3)

		# Get a fresh copy of test and train graphs for the new epoch.
		graphs_pos_train , graphs_pos_train_epoch = tee(graphs_pos_train)
		graphs_neg_train , graphs_neg_train_epoch = tee(graphs_neg_train)

		# Take n% of test and train graph sets. n is the epoch counter.
		graphs_pos_train_epoch , iterable_not_used = random_bipartition_iter(graphs_pos_train_epoch , relative_size=relative_size)

		graphs_neg_train_epoch , iterable_not_used = random_bipartition_iter(graphs_neg_train_epoch , relative_size=relative_size)

		# Make copies of each for the whole epoch run.
		graphs_pos_train_epoch_1 , graphs_pos_train_epoch_2, graphs_pos_train_epoch_3  = tee(graphs_pos_train_epoch, 3)

		graphs_neg_train_epoch_1 , graphs_neg_train_epoch_2 = tee(graphs_neg_train_epoch)
		
		# Instantiate the vectorizer.
		vectorizer = Vectorizer( complexity = 2 )

		# Train TrueSamplesModel classifier.
		estimator = fit(graphs_pos_train_epoch_1, graphs_neg_train_epoch_1, vectorizer, n_jobs=1, cv=3)

		# Design the batch of new sequences.
		designed_iterable = produce_batch_RNA(graphs = graphs_pos_train_epoch_2, vectorizer = vectorizer, estimator = estimator, **opts)

		# Make copies of designed_iterable.
		designed_iterable , designed_iterable2 =  tee(designed_iterable)
		
		# Build the graphs iterator over the designed sequences.
		graphs_synthesized = rnafold_to_eden( designed_iterable )

		# Produce graph iterator over negative samples for synthesized data.
		iterable_neg = seq_to_seq( designed_iterable2, modifier=shuffle_modifier, times=3, order=2 )
		graphs_neg_synthesized = rnafold_to_eden( iterable_neg )

		# Mix the sample with "true sequences".
		graphs_mixed_pos = chain(graphs_pos_train_epoch_3, graphs_synthesized)
		graphs_mixed_neg = chain(graphs_neg_train_epoch_2, graphs_neg_synthesized)

		# Train MixedSamplesModel classifier.
		estimator2 = fit(graphs_mixed_pos, graphs_mixed_neg, vectorizer, n_jobs=1, cv=3)

		# Test the test set against TrueSamplesModel -> output performance
		logger.info('Evaluating the True model performance:')
		roc_t , apr_t = estimate(graphs_pos_test_epoch_1, graphs_neg_test_epoch_1 , estimator, vectorizer, n_jobs=1)

		# Test the test set against SynthesizedSamplesModel -> output performance
		logger.info('Evaluating the mixed model performance:')
		roc_s , apr_s = estimate(graphs_pos_test_epoch_2, graphs_neg_test_epoch_2, estimator2, vectorizer, n_jobs=1)
		
		## Additional sample counts:
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
	parser.add_argument('--rfam_id', '-i', type = str , default = 'RF00005')
	parser.add_argument('--log_file', '-l', type = str , default = '~/Synthesis.log')
	parser.add_argument('--antaRNA_params', '-p', type = str , default = './antaRNA.ini')
	parser.add_argument('--nt_importance_threshold', '-a', type = int , default = 0)
	parser.add_argument('--nmin_important_nt_adjaceny', '-b', type = int , default = 1)
	parser.add_argument('--bp_importance_threshold', '-c', type = int , default = 0)
	parser.add_argument('--nmin_important_bp_adjaceny', '-d', type = int , default = 1)
	parser.add_argument('--nmin_unpaired_nt_adjacency', '-e', type = int , default = 1)
	parser.add_argument('--multi_sequence_size', '-n', type = int , default = 1)
	parser.add_argument('--filtering_threshold', '-f', type = int , default = 0)
	parser.add_argument('--batch_proportion', '-g', type = int , default = 1)
	parser.add_argument('--epoch_instances', '-k', type = int , default = 10)
	parser.add_argument('--experiment_runs', '-j', type = int , default = 10)
	parser.add_argument('--split_ratio', '-r' , type = float , default = 0.2)
	args = parser.parse_args()

	params = {'rfam_id':args.rfam_id , 'log_file':args.log_file , 'antaRNA_param_file':args.antaRNA_params , \
			'nt_importance_threshold':args.nt_importance_threshold , 'nmin_important_nt_adjaceny':args.nmin_important_nt_adjaceny ,\
			'bp_importance_threshold':args.bp_importance_threshold , 'nmin_important_bp_adjaceny':args.nmin_important_bp_adjaceny ,\
			'nmin_unpaired_nt_adjacency':args.nmin_unpaired_nt_adjacency , 'multi_sequence_size':args.multi_sequence_size ,\
			'filtering_threshold':args.filtering_threshold , 'batch_proportion':args.batch_proportion , 'epoch_instances':args.epoch_instances ,\
			'experiment_runs':args.experiment_runs, 'split_ratio':args.split_ratio}

	return params


def experiment():
	"""
	Main body of RNASynthesis experiment.
	"""
	#logger = logging.getLogger()
	#logger.basicConfig(level=logging.INFO)

	params = get_args()

	iterable_pos = fasta_to_sequence( rfam_url(params['rfam_id']) )
	iterable_pos, iterable_pos_ = tee(iterable_pos)
	iterable_neg = seq_to_seq( iterable_pos_ , modifier=shuffle_modifier , times=3 , order=2 )
	
	#Positive sample graphs.
	graphs_pos = rnafold_to_eden(iterable_pos)
	#Negative sample graphs.
	graphs_neg = rnafold_to_eden(iterable_neg)
	#Split samples to train and test iterables.
	graphs_pos_train, graphs_pos_test = random_bipartition_iter(graphs_pos, relative_size= params['split_ratio'])
	graphs_neg_train, graphs_neg_test = random_bipartition_iter(graphs_neg, relative_size= params['split_ratio'])

	exp_roc_t = []
	exp_roc_s = []
	exp_apr_t = []
	exp_apr_s = []

	for counter in range(1 , 11):
		logger.info('Starting epoch %s:' %counter)
		graphs_pos_train,graphs_pos_train_c = tee(graphs_pos_train)
		graphs_pos_test,graphs_pos_test_c = tee(graphs_pos_test)
		graphs_neg_train,graphs_neg_train_c = tee(graphs_neg_train)
		graphs_neg_test,graphs_neg_test_c = tee(graphs_neg_test)
		
		prct = float(counter)/10
		
		mroct, maprt, mrocs, maprs, elapsed_time = run_epoch(graphs_pos_test =graphs_pos_test_c , graphs_neg_test = graphs_neg_test_c ,\
							graphs_pos_train = graphs_pos_train_c , graphs_neg_train = graphs_neg_train_c , \
							nt_importance_threshold= params['nt_importance_threshold'], nmin_important_nt_adjaceny=params['nmin_important_nt_adjaceny'], \
							bp_importance_threshold=params['bp_importance_threshold'], nmin_important_bp_adjaceny=params['nmin_important_bp_adjaceny'], \
							nmin_unpaired_nt_adjacency=params['nmin_unpaired_nt_adjacency'], synthesized_batch_proportion=params['batch_proportion'], \
							multi_sequence_size=params['multi_sequence_size'], filtering_threshold =params['filtering_threshold'], \
							runs = params['epoch_instances'], relative_size = prct, antaRNA_param_file =params['antaRNA_param_file'] ,\
							performance_log_file =params['log_file'])

		sys.stdout.write('Performance measures for epoch: %d\n' %counter)
		sys.stdout.write('ROC for True samples: \n')
		sys.stdout.write('%s \n' %mroct)
		sys.stdout.write('ROC for Mixed samples: \n')
		sys.stdout.write('%s \n' %mrocs)
		sys.stdout.write('APR for True samples: \n')
		sys.stdout.write('%s \n' %maprt)
		sys.stdout.write('APR for Mixed samples: \n')
		sys.stdout.write('%s \n' %maprs)
		exp_roc_t.append(mroct)
		exp_roc_s.append(mrocs)
		exp_apr_t.append(maprt)
		exp_apr_s.append(maprs)
	return 	exp_roc_t.append(mroct) , exp_roc_s.append(mrocs) , exp_apr_t.append(maprt) , exp_apr_s.append(maprs)



if __name__ == "__main__":
	"""
	Main body of the experiment.
	Runs on a specific rfam family.
	Starts from a portion of the training samples and increases the portion per iteration.
	Runs an epoch on each portion.
	Saves the results in individual files labeld with the portion.
	"""
	logging.basicConfig(level=logging.INFO)
	logger = logging.getLogger(__name__)
	logger.info('Call to experiment module.')


	experiment()

