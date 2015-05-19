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


def output_measure_column(measure_list, measure , header = ''):
	"""
	Output formatting function.
	"""
	output = ''
	measure_column = ''

	if not header == '':
		output = header + '\n'
	output = output + measure + '\n'

	for measure in measure_list:
		measure_column = measure_column + str(measure) + '\n'
		
	return output + measure_column


def flush_to_file(iterable, file_name, header =''):
	"""
	Generic function
	Outputs test results into file.
	"""
	with open(file_name,'a') as f:
		if not header=='':
			f.write(header)
			f.write('\n')
		for item in iterable:
			f.write(item)
			f.write('\n')
	return


def output_results(opts=None, mroct=None, mrocs=None, maprt=None, maprs=None):
	"""
	Outputs test results into file.
	This function is written specifically for this experiment.
	"""
	performance_log_file = opts.get('performance_log_file')
	relative_size = opts.get('relative_size')
	header = 'RNADesign test on rfam family %s. \n' %rfam_id
	header = header + 'Performance measure collection on 10 runs of %f of data.\n' %relative_size
	header = header + 'Overall time %f seconds. \n' %elapsed_time

	header = header + 'Sequence design constraints for this experiment: \n'
	Cstr_threshold = opts.get('nt_importance_threshold')
	Cstr_adjaceny = opts.get('nmin_important_nt_adjaceny')
	DotBr_threshold = opts.get('bp_importance_threshold')
	DotBr_adjacency = opts.get('nmin_important_bp_adjaceny')
	unpaired_adjacency = opts.get('nmin_unpaired_nt_adjacency')

	header = header + 'Generic model sequence constraints: \n'
	header = header + 'Nucleotide importance threshold:\t\t%f\n' %Cstr_threshold
	header = header + 'Minimum number of adjacent important nucleotides:\t\t%d\n' %Cstr_adjaceny

	header = header + 'Dot-bracket notation constraints: \n'
	header = header + 'Nucleotide importance threshold:\t\t%f\n' %DotBr_threshold
	header = header + 'Minimum number of adjacent base pairs:\t\t%d\n' %DotBr_adjacency

	header = header + 'Unpaired regions constraint: \n'
	header = header + 'Minimum number of adjacent unpaired nodes:\t\t%d\n' %unpaired_adjacency

	measures_list = []

	measures_list.append( output_measure_column(mroct, 'ROC' , header = 'True Samples' ))
	measures_list.append( output_measure_column(mrocs, 'ROC' , header = 'Mixed Samples' ))
	measures_list.append( output_measure_column(maprt, 'APR' , header = 'True Samples' ))
	measures_list.append( output_measure_column(maprs, 'APR' , header = 'Mixed Samples' ))

	flush_to_file(measures_list , performance_log_file, header)


def rfam_url(family_id):
	return 'http://rfam.xfam.org/family/%s/alignment?acc=%s&format=fastau&download=0'%(family_id,family_id)


def run_epoch(runs = None , antaRNA_param_file = None , performance_log_file = None , graphs_pos_test = None , graphs_neg_test = None , \
		graphs_pos_train = None , graphs_neg_train = None , nt_importance_threshold=0, nmin_important_nt_adjaceny=1, \
		bp_importance_threshold=0, nmin_important_bp_adjaceny=1, nmin_unpaired_nt_adjacency=1, \
		synthesized_batch_proportion=1, multi_sequence_size=1, filtering_threshold = 0 , relative_size = None, rfam_id=None):
	"""
	Executes one epoch of n runs.
	Each run does the experiment on "relative_size" of training samples.
	Outputs four lists containing the overall ROC and APR performance measures for true and mixed training samples.
	Also outputs the overall elapsed time on n runs.
	"""
	start_time = time.time()
	
	measure_list_ROCT = []
	measure_list_APRT = []
	measure_list_ROCS = []
	measure_list_APRS = []
	
	for epoch in range(1,runs + 1):
		
		sequence_pool = create_sequence_pool(rfam_id)
		sys.stdout.write('Initial lengt of the sequence pool: %d\n' %len(sequence_pool)) 
				
		graphs_pos_test, graphs_pos_test_epoch_1, graphs_pos_test_epoch_2 = tee(graphs_pos_test,3)
		graphs_neg_test, graphs_neg_test_epoch_1, graphs_neg_test_epoch_2 = tee(graphs_neg_test,3)

		# Get a fresh copy of test and train graphs for the new epoch.
		graphs_pos_train , graphs_pos_train_epoch = tee(graphs_pos_train)
		graphs_neg_train , graphs_neg_train_epoch = tee(graphs_neg_train)

		# Take n% of test and train graph sets. n is the epoch counter.
		graphs_pos_train_epoch , iterable_not_used = random_bipartition_iter(graphs_pos_train_epoch , relative_size=relative_size)

		graphs_neg_train_epoch , iterable_not_used = random_bipartition_iter(graphs_neg_train_epoch , relative_size=relative_size)

		# Make m copies of each for the whole epoch run.
		graphs_pos_train_epoch_1 , graphs_pos_train_epoch_2, graphs_pos_train_epoch_3, graphs_pos_train_epoch_counter  = tee(graphs_pos_train_epoch, 4)

		graphs_neg_train_epoch_1 , graphs_neg_train_epoch_2 = tee(graphs_neg_train_epoch)
		
		# Calculate the batch size proportional to the size of the epoch's training samples.
		epoch_train_pos_size = 0
		for graph in graphs_pos_train_epoch_counter:
			epoch_train_pos_size += 1
		# Instantiate the vectorizer.
		vectorizer = Vectorizer( complexity=2 )

		# Train TrueSamplesModel classifier.
		estimator = fit(graphs_pos_train_epoch_1, graphs_neg_train_epoch_1, vectorizer, n_jobs=1, cv=3)


		sys.stdout.write('Number of training samples: %d\n' %epoch_train_pos_size)
		batch_size = synthesized_batch_proportion * epoch_train_pos_size
		sys.stdout.write('Number of samples to be synthesized:: %d\n' %batch_size)
			
		# Design the batch of new sequences.
		designed_sequence_list, time_elapsed = rdt.design_batch_RNA(param_file=antaRNA_param_file , iterable=graphs_pos_train_epoch_2, vectorizer=vectorizer,\
												design_estimator=estimator, filter_estimator=estimator , nt_importance_threshold=nt_importance_threshold, \
												nmin_important_nt_adjaceny=nmin_important_nt_adjaceny, bp_importance_threshold=bp_importance_threshold, \
												nmin_important_bp_adjaceny=nmin_important_bp_adjaceny, nmin_unpaired_nt_adjacency=nmin_unpaired_nt_adjacency,\
												batch_size=batch_size, multi_sequence_size=multi_sequence_size, filtering_threshold = filtering_threshold,\
												sequence_pool=sequence_pool)

		#sys.stdout.write('Lengt of the sequence pool after batch design: %d\n' %len(sequence_pool))

		if len(designed_sequence_list) <= 0:
			logger.info('The application was unable to produce new sequences in epoch %d. Quitting this run...' %epoch)
			break
		else:
			logger.info('The application produced %d new sequences in epoch %d.' %(len(designed_sequence_list)/2, epoch))
			
		# Build the graph list of designed sequences.
		iterable_pos = fasta_to_sequence( designed_sequence_list )
		graphs_synthesized = rnafold_to_eden( iterable_pos )

		iterable_neg = fasta_to_sequence( designed_sequence_list, modifier=shuffle_modifier, times=3, order=2 )
		graphs_neg_synthesized = rnafold_to_eden( iterable_pos )

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
		
		measure_list_ROCT.append((roc_t,epoch_train_pos_size,epoch_train_pos_size))
		measure_list_APRT.append((apr_t,epoch_train_pos_size,epoch_train_pos_size))
		measure_list_ROCS.append((roc_s, batch_size + epoch_train_pos_size , len(designed_sequence_list)/2 + epoch_train_pos_size))
		measure_list_APRS.append((apr_s, batch_size + epoch_train_pos_size , len(designed_sequence_list)/2 + epoch_train_pos_size))
	
	elapsed_time = time.time() - start_time
	
	return measure_list_ROCT, measure_list_APRT, measure_list_ROCS, measure_list_APRS, elapsed_time


if __name__ == "__main__":
	"""
	Main body of the experiment.
	Runs on a specific rfam family.
	Starts from a portion of the training samples and increases the portion per iteration.
	Runs an epoch on each portion.
	Saves the results in individual files labeld with the portion.
	"""
	start_portion = 1
	end_portion = 11
	
	logging.basicConfig(level=logging.INFO)
	logger = logging.getLogger(__name__)
	logger.info('Call to experiment module.')
	
	rfam_id = 'RF00005'
	opts={'runs':2 , 'antaRNA_param_file':'/home/kohvaeip/RLS/lib/antaRNA/antaRNA.ini' , \
	  'performance_log_file':'/home/kohvaeip/perf.log', 'nt_importance_threshold':-1.1 , 'nmin_important_nt_adjaceny':1 , \
	  'bp_importance_threshold':-0.85 , 'nmin_important_bp_adjaceny':3 , 'nmin_unpaired_nt_adjacency':3 , \
	  'synthesized_batch_proportion':2, 'multi_sequence_size':3, 'filtering_threshold':0, 'relative_size':0.3, 'rfam_id':rfam_id}
	  
	iterable_pos = fasta_to_sequence( rfam_url(rfam_id) )
	iterable_pos, iterable_pos_ = tee(iterable_pos)
	iterable_neg = seq_to_seq( iterable_pos_ , modifier=shuffle_modifier , times=3 , order=2 )

	#Positive sample graphs.
	graphs_pos = rnafold_to_eden(iterable_pos)
	#Negative sample graphs.
	graphs_neg = rnafold_to_eden(iterable_neg)
	#Keep 0.8 of data samples for testing.
	graphs_pos_train, graphs_pos_test = random_bipartition_iter(graphs_pos, relative_size=0.2)
	graphs_neg_train, graphs_neg_test = random_bipartition_iter(graphs_neg, relative_size=0.2)
	
	for counter in range(start_portion , end_portion):
		logger.info('Starting batch run number %s:' %counter)
		graphs_pos_train,graphs_pos_train_c = tee(graphs_pos_train)
		graphs_pos_test,graphs_pos_test_c = tee(graphs_pos_test)
		graphs_neg_train,graphs_neg_train_c = tee(graphs_neg_train)
		graphs_neg_test,graphs_neg_test_c = tee(graphs_neg_test)
		
		prct = float(counter)/10
		opts['performance_log_file']='/home/kohvaeip/perf.log'+str(prct)
		opts['relative_size']=prct
		mroct, maprt, mrocs, maprs, elapsed_time = run_epoch(graphs_pos_test =graphs_pos_test_c , graphs_neg_test = graphs_neg_test_c ,\
							graphs_pos_train = graphs_pos_train_c , graphs_neg_train = graphs_neg_train_c , **opts)
		if len(mroct) > 0:
			output_results(opts=opts, mroct=mroct, mrocs=mrocs, maprt=maprt, maprs=maprs)
