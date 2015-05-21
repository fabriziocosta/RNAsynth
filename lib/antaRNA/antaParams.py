#!/usr/bin/env python

import os
import sys
import ConfigParser
import argparse
import logging

'''
structure = ''
Cseq = ''
tGC = 0.5 removed
colonies = 10
name = "antaRNA_"
alpha = 1.0
beta = 1.0
evaporation_rate = 0.2
struct_correction_term = 0.5
GC_correction_term = 5.0
seq_correction_term = 1.0
degreeOfSequenceInducement = 1
file_id = "STDOUT"
verbose = ''
output_verbose = True
tGCmax = -1.0
tGCvar = -1.0
termination_convergence = 50
convergence_count = 130
reset_limit = 5
improve = "s"
seed = "none"
temperature = 37.0
paramFile = ""
return_mod = True
'''

defaults = {'structure':'', 'Cseq':'','colonies':'10', 'name':'antaRNA_', 'alpha':'1.0', \
			'beta':'1.0', 'evaporation_rate':'0.2', 'struct_correction_term':'0.5', 'GC_correction_term':'5.0', \
			'seq_correction_term':'1.0', 'degreeOfSequenceInducement':'1', 'file_id':'STDOUT', 'verbose':'', \
			'output_verbose':'True', 'tGCmax':'-1.0', 'tGCvar':'-1.0', 'termination_convergence':'50', 'convergence_count':'130', \
			'reset_limit':'5', 'improve':'s', 'seed':'none', 'temperature':'37.0', 'paramFile':'', 'return_mod':'True'}


class antaParams():
	
	def __init__(self, param_file):
		try:
			config = ConfigParser.SafeConfigParser(defaults)
			config.read(param_file)
		   
			#self.tGC = float(config.get('PARAMETERS', 'tGC')) if config.get('PARAMETERS', 'tGC') else 0.5
			self.colonies = int(config.get('PARAMETERS', 'colonies')) if config.get('PARAMETERS', 'colonies') else 10
			self.name = config.get('PARAMETERS', 'name') if config.get('PARAMETERS', 'name') else 'antaRNA_'
			self.alpha = float(config.get('PARAMETERS', 'alpha')) if config.get('PARAMETERS', 'alpha') else 1.0
			self.beta = float(config.get('PARAMETERS', 'beta')) if config.get('PARAMETERS', 'beta') else 1.0
			self.evaporation_rate = float(config.get('PARAMETERS', 'evaporation_rate')) if config.get('PARAMETERS', 'evaporation_rate') else 0.2
			self.struct_correction_term = float(config.get('PARAMETERS', 'struct_correction_term')) if config.get('PARAMETERS', 'struct_correction_term') else 0.5
			self.GC_correction_term = float(config.get('PARAMETERS', 'GC_correction_term')) if config.get('PARAMETERS', 'GC_correction_term') else 5.0
			self.seq_correction_term = float(config.get('PARAMETERS', 'seq_correction_term')) if config.get('PARAMETERS', 'seq_correction_term') else 1.0
			self.degreeOfSequenceInducement = int(config.get('PARAMETERS', 'degreeOfSequenceInducement')) if config.get('PARAMETERS', 'degreeOfSequenceInducement') else 1
			self.file_id = config.get('PARAMETERS', 'file_id') if config.get('PARAMETERS', 'file_id') else 'STDOUT'
			self.verbose = bool(config.get('PARAMETERS', 'verbose')) if config.get('PARAMETERS', 'verbose') else False
			self.output_verbose = bool(config.get('PARAMETERS', 'output_verbose')) if config.get('PARAMETERS', 'output_verbose') else False
			self.tGCmax = float(config.get('PARAMETERS', 'tGCmax')) if config.get('PARAMETERS', 'tGCmax') else -1.0
			self.tGCvar = float(config.get('PARAMETERS', 'tGCvar')) if config.get('PARAMETERS', 'tGCvar') else -1.0
			self.termination_convergence = int(config.get('PARAMETERS', 'termination_convergence')) if config.get('PARAMETERS', 'termination_convergence') else 50
			self.convergence_count = int(config.get('PARAMETERS', 'convergence_count')) if config.get('PARAMETERS', 'convergence_count') else 130
			self.reset_limit = int(config.get('PARAMETERS', 'reset_limit')) if config.get('PARAMETERS', 'reset_limit') else 5
			self.improve = config.get('PARAMETERS', 'improve') if config.get('PARAMETERS', 'improve') else 's'
			self.temperature = float(config.get('PARAMETERS', 'temperature')) if config.get('PARAMETERS', 'temperature') else 37.0
			self.paramFile = config.get('PARAMETERS', 'paramFile') if config.get('PARAMETERS', 'paramFile') else ''
			self.return_mod = bool(config.get('PARAMETERS', 'return_mod')) if config.get('PARAMETERS', 'return_mod') else True
			self.seed = config.get('PARAMETERS', 'seed') if config.get('PARAMETERS', 'seed') else 'none'

		except Exception as e: 
			print e
		#except:
			#logger.error('Failed to initialize parameter container.') 
			#sys.exit()
			
	def to_dict(self):
		return {'paramFile': self.paramFile, 'colonies': self.colonies, 'name': self.name, 'alpha': self.alpha, \
				'beta': self.beta, 'evaporation_rate': self.evaporation_rate, 'struct_correction_term': self.struct_correction_term, \
				'GC_correction_term': self.GC_correction_term, 'seq_correction_term': self.seq_correction_term, \
				'degreeOfSequenceInducement': self.degreeOfSequenceInducement, 'file_id': self.file_id, 'verbose': self.verbose, \
				'output_verbose': self.output_verbose, 'tGCmax': self.tGCmax, 'tGCvar': self.tGCvar, 'termination_convergence': self.termination_convergence, \
				'convergence_count': self.convergence_count, 'reset_limit': self.reset_limit, 'improve': self.improve, 'seed': self.seed, \
				'temperature': self.temperature, 'return_mod': self.return_mod}


if __name__ == "__main__":
		
	logging.basicConfig(level=logging.INFO)
	logger = logging.getLogger(__name__)
	logger.info('antaRNA parameter container called.')

 

