#!/usr/bin/env python

import logging
import sys

from lib.antaRNA import antaRNA_v117
from lib.antaRNA import antaRNA_v109

logger = logging.getLogger(__name__)


class AbstractDesigner(object):

    """Interface declaration for designer classes."""

    def design(self, constraints=None):
        raise NotImplementedError("Design method not implemented.")


class AntaRNAv109Designer(AbstractDesigner):

    def __init__(self,
                 tGC=0.5,
                 colonies=1,
                 name='antaRNA_',
                 alpha=1.0,
                 beta=1.0,
                 evaporation_rate=0.2,
                 struct_correction_term=0.5,
                 GC_correction_term=5.0,
                 seq_correction_term=1.0,
                 degreeOfSequenceInducement=1,
                 file_id='STDOUT',
                 verbose=False,
                 output_verbose=False,
                 tGCmax=-1.0,
                 tGCvar=-1.0,
                 termination_convergence=50,
                 convergence_count=130,
                 reset_limit=5,
                 improve='s',
                 temperature=37.0,
                 paramFile='',
                 return_mod=True,
                 seed='none'
                 ):
        self.params = {}
        self.params['tGC'] = tGC
        self.params['colonies'] = colonies
        self.params['name'] = 'antaRNA_'
        self.params['alpha'] = alpha
        self.params['beta'] = beta
        self.params['evaporation_rate'] = evaporation_rate
        self.params['struct_correction_term'] = struct_correction_term
        self.params['GC_correction_term'] = GC_correction_term
        self.params['seq_correction_term'] = seq_correction_term
        self.params['degreeOfSequenceInducement'] = degreeOfSequenceInducement
        self.params['file_id'] = file_id
        self.params['verbose'] = verbose
        self.params['output_verbose'] = output_verbose
        self.params['tGCmax'] = tGCmax
        self.params['tGCvar'] = tGCvar
        self.params['termination_convergence'] = termination_convergence
        self.params['convergence_count'] = convergence_count
        self.params['reset_limit'] = reset_limit
        self.params['improve'] = improve
        self.params['temperature'] = temperature
        self.params['paramFile'] = paramFile
        self.params['return_mod'] = return_mod
        self.params['seed'] = seed
        logger.info('Instantiated an instance of AntaRNAv109Designer.')

    # dot_bracket_constraint_string, sequence_constraint_string):
    def design(self, constraints=None):
        dot_bracket_constraint_string = constraints[0]
        sequence_constraint_string = constraints[1]
        self.params['tGC'] = constraints[2]
        result = ', '.join(antaRNA_v109.findSequence(
            dot_bracket_constraint_string, sequence_constraint_string, **self.params))
        sequence = result.split("\n")[2]
        return sequence


class AntaRNAv117Designer(AbstractDesigner):

    def __init__(self,
                 Cstr="",
                 Cseq="",
                 tGC=[],
                 level=1,
                 tGCmax=-1.0,
                 tGCvar=-1.0,
                 temperature=37.0,
                 paramFile="",
                 noGUBasePair=False,
                 noLBPmanagement=True,
                 pseudoknots=False,
                 usedProgram="RNAfold",
                 pkprogram="pKiss",
                 pkparameter=False,
                 HotKnots_PATH="",
                 strategy="A",
                 noOfColonies=1,
                 output_file="STDOUT",
                 py=True,
                 name="antaRNA",
                 verbose=False,
                 output_verbose=False,
                 seed="none",
                 improve_procedure="s",
                 Resets=5,
                 ants_per_selection=10,
                 ConvergenceCount=130,
                 antsTerConv=50,
                 alpha=1.0,
                 beta=1.0,
                 ER=0.2,
                 Cstrweight=0.5,
                 Cgcweight=5.0,
                 Cseqweight=1.0,
                 omega=2.23,
                 time=600):

        self.designer = antaRNA_v117.AntHill()
        self.designer.params.Cstr = Cstr
        self.designer.params.Cseq = Cseq
        self.designer.params.tGC = tGC
        self.designer.params.level = level
        self.designer.params.tGCmax = tGCmax
        self.designer.params.tGCvar = tGCvar
        self.designer.params.temperature = temperature
        self.designer.params.paramFile = paramFile
        self.designer.params.noGUBasePair = noGUBasePair
        self.designer.params.noLBPmanagement = noLBPmanagement
        self.designer.params.pseudoknots = pseudoknots
        self.designer.params.pkprogram = pkprogram
        self.designer.params.pkparameter = pkparameter
        self.designer.params.HotKnots_PATH = HotKnots_PATH
        self.designer.params.strategy = strategy
        self.designer.params.noOfColonies = noOfColonies
        self.designer.params.output_file = output_file
        self.designer.params.py = py
        self.designer.params.name = name
        self.designer.params.verbose = verbose
        self.designer.params.output_verbose = output_verbose
        self.designer.params.seed = seed
        self.designer.params.improve_procedure = improve_procedure
        self.designer.params.Resets = Resets
        self.designer.params.ants_per_selection = ants_per_selection
        self.designer.params.ConvergenceCount = ConvergenceCount
        self.designer.params.antsTerConv = antsTerConv
        self.designer.params.alpha = alpha
        self.designer.params.beta = beta
        self.designer.params.ER = ER
        self.designer.params.Cstrweight = Cstrweight
        self.designer.params.Cgcweight = Cgcweight
        self.designer.params.Cseqweight = Cseqweight
        self.designer.params.omega = omega
        self.designer.params.time = time

        if not(self.designer.params.error == "0"):
            logger.info(
                'antaRNA parameters integrity check faild with the following error:')
            logger.info('%s' % self.designer.params.error)
            sys.exit()

        logger.info('Instantiated an instance of AntaRNAv117Designer.')

    def design(self, constraints=None):
        """
        DOCUMENTATION
        """
        # check first
        self.designer.params.Cstr = constraints[0]
        self.designer.params.Cseq = constraints[1]
        if type(constraints[2]) == list:
            self.designer.params.tGC = constraints[2]
        else:
            self.designer.params.tGC = [constraints[2]]
        self.designer.params.check()
        self.designer.swarm()
        result = ''
        for r in self.designer.result:
            result = r[1].split(":")[1]
        return result


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to rna_desiner module.')
    constraints = ("...(((...)))...((...))...",
                   "NNNNNNNNNNNNNNNNNNNNNNNNU", [0.5])
    designer = AntaRNAv117Designer()

    r = designer.design(constraints=constraints)
    print r

    # designer2 = AntaRNAv109Designer()
    # result = designer2.design(constraints=("...(((...)))...((...))...",
    #                                        "NNNNNNNNNNNNNNNNNNNNNNNNU", 0.3))
