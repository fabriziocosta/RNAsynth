#!/usr/bin/env python

import logging
import argparse
import sys

from lib.antaRNA import antaRNA_v109

logger = logging.getLogger(__name__)


class AbstractDesigner(object):

    """Interface declaration for designer classes."""

    def design(self, constraints=None):
        raise NotImplementedError("Design method not implemented.")


class antRNADesigner(AbstractDesigner):

    def __init__(self,
				 Cstr = "",
				 Cseq = "",
				 tGC = [],
				 level = 1,
				 tGCmax = -1.0,
				 tGCvar = -1.0,
				 temperature = 37.0,
				 paramFile = "",
				 noGUBasePair = False,
				 noLBPmanagement = True,
				 pseudoknots = False,
				 usedProgram = "RNAfold",
				 pkprogram = "pKiss",
				 pkparameter = False,
				 HotKnots_PATH = "",
				 strategy = "A",
				 noOfColonies = 1,
				 output_file = "STDOUT",
				 py = True,
				 name="antaRNA",
				 verbose = False ,
				 output_verbose = False,
				 seed ="none",
				 improve_procedure = "s",
				 Resets = 5,
				 ants_per_selection = 10,
				 ConvergenceCount = 130,
				 antsTerConv = 50,
				 alpha = 1.0,
				 beta = 1.0,
				 ER = 0.2,
				 Cstrweight = 0.5,
				 Cgcweight = 5.0,
				 Cseqweight = 1.0,
				 omega = 2.23,
				 time = 600)
		
        self.designer = AntHill()
        self.designer.Cstr = Cstr
		self.designer.Cseq = Cseq
		self.designer.tGC = tGC
		self.designer.level = level
		self.designer.tGCmax = tGCmax
		self.designer.tGCvar = tGCvar
		self.designer.temperature = temperature
		self.designer.paramFile = paramFile
		self.designer.noGUBasePair = noGUBasePair
		self.designer.noLBPmanagement = noLBPmanagement
		self.designer.pseudoknots = pseudoknots
		self.designer.pkprogram = pkprogram
		self.designer.pkparameter = pkparameter
		self.designer.HotKnots_PATH = HotKnots_PATH
		self.designer.strategy = strategy
		self.designer.noOfColonies = noOfColonies
		self.designer.output_file = output_file
		self.designer.py = py
		self.designer.name = name
		self.designer.verbose = verbose
		self.designer.output_verbose = output_verbose
		self.designer.seed = seed
		self.designer.improve_procedure = improve_procedure
		self.designer.Resets = Resets
		self.ants_per_selection = ants_per_selection
		self.designer.ConvergenceCount = ConvergenceCount
		self.designer.antsTerConv = antsTerConv
		self.designer.alpha = alpha
		self.designer.beta = beta
		self.designer.ER = ER
		self.designer.Cstrweight = Cstrweight
		self.designer.Cgcweight = Cgcweight
		self.designer.Cseqweight = Cseqweight
		self.designer.omega = omega
		self.designer.time = time
        # self.designer.params.py = False
        self.designer.params.check()
        if not(self.designer.params.error == "0"):
			logger.info('antaRNA parameters integrity check faild with the following error:')
			logger.info('%s' %self.designer.params.error)
			sys.exit()
			
        logger.info('Instantiated an instance of antRNADesigner object.')

    def design(self, constraints=None):
        """
        DOCUMENTATION
        """
        # check first
        # self.designer.Cstr = constraints[0]
        # self.designer.Cseq = constraints[1]
		# self.designer.tGC = constraints[2]
        result = self.designer.swarm()
		# extract generated sequence from output. --> sequence
        return sequence


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to rna_desiner module.')
    designer = RNADesign()
    seq = designer.design(
        "...(((...)))...((...))...", "NNNNNNNNNNNNNNNNNNNNNNNNU")
    print seq
