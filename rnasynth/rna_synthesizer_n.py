#!/usr/bin/env python


from itertools import tee
from itertools import izip
import logging

from sklearn.linear_model import SGDClassifier

from eden.converter.fasta import fasta_to_sequence
from eden.converter.rna.rnashapes import rnashapes_to_eden
from eden.util import fit
from eden.graph import Vectorizer
from eden.modifier.seq import seq_to_seq
from eden.modifier.seq import shuffle_modifier

from rnasynth.constraint_extractor import ConstraintExtractor
from rnasynth.rna_designer import AntaRNAv117Designer


logger = logging.getLogger(__name__)


class Initializer():

    def __init__(self,
                 # constraint_extractor params
                 importance_threshold_sequence_constraint=0,
                 min_size_connected_component_sequence_constraint=1,
                 importance_threshold_structure_constraint=0,
                 min_size_connected_component_structure_constraint=1,
                 min_size_connected_component_unpaired_structure_constraint=1,
                 # antarna_designer_v117 params
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
                 time=600,
                 # pre_process params
                 shape_type=5,
                 energy_range=35,
                 max_num=3,
                 split_components=True,
                 # rna_synthesizer params
                 instance_score_threshold=0,
                 shuffle_order=2,
                 negative_shuffle_ratio=2,
                 vectorizer_complexity=2,
                 max_num_graphs_per_seq=3,
                 n_jobs=-1,
                 cv=3,
                 n_iter_search=1,
                 n_synthesized_sequences_per_seed_sequence=3
                 ):

        self.constraint_extractor = ConstraintExtractor(
            importance_threshold_sequence_constraint=importance_threshold_sequence_constraint,
            min_size_connected_component_sequence_constraint=min_size_connected_component_sequence_constraint,
            importance_threshold_structure_constraint=importance_threshold_structure_constraint,
            min_size_connected_component_structure_constraint=min_size_connected_component_structure_constraint,
            min_size_connected_component_unpaired_structure_constraint=min_size_connected_component_unpaired_structure_constraint)

        self.designer = AntaRNAv117Designer(Cstr=Cstr,
                                            Cseq=Cseq,
                                            tGC=tGC,
                                            level=level,
                                            tGCmax=tGCmax,
                                            tGCvar=tGCvar,
                                            temperature=temperature,
                                            paramFile=paramFile,
                                            noGUBasePair=noGUBasePair,
                                            noLBPmanagement=noLBPmanagement,
                                            pseudoknots=pseudoknots,
                                            usedProgram=usedProgram,
                                            pkprogram=pkprogram,
                                            pkparameter=pkparameter,
                                            HotKnots_PATH=HotKnots_PATH,
                                            strategy=strategy,
                                            noOfColonies=noOfColonies,
                                            output_file=output_file,
                                            py=py,
                                            name=name,
                                            verbose=verbose,
                                            output_verbose=output_verbose,
                                            seed=seed,
                                            improve_procedure=improve_procedure,
                                            Resets=Resets,
                                            ants_per_selection=ants_per_selection,
                                            ConvergenceCount=ConvergenceCount,
                                            antsTerConv=antsTerConv,
                                            alpha=alpha,
                                            beta=beta,
                                            ER=ER,
                                            Cstrweight=Cstrweight,
                                            Cgcweight=Cgcweight,
                                            Cseqweight=Cseqweight,
                                            omega=omega,
                                            time=time)

        self.pre_processor = PreProcessor(shape_type=shape_type,
                                          energy_range=energy_range,
                                          max_num=max_num,
                                          split_components=split_components)

        self.estimator = SGDClassifier()
        self.vectorizer = Vectorizer()
        self.synthesizer = RNASynth(estimator=self.estimator,
                                    vectorizer=self.vectorizer,
                                    designer=self.designer,
                                    pre_processor=self.pre_processor,
                                    constraint_extractor=self.constraint_extractor,
                                    n_synthesized_sequences_per_seed_sequence=n_synthesized_sequences_per_seed_sequence,
                                    instance_score_threshold=instance_score_threshold,
                                    shuffle_order=shuffle_order,
                                    negative_shuffle_ratio=negative_shuffle_ratio,
                                    n_jobs=n_jobs,
                                    cv=cv,
                                    n_iter_search=n_iter_search)
        logger.info('Created a RNASynthesizer object.')

    def init_synthesizer(self):
        return self.synthesizer


class PreProcessor():

    def __init__(self,
                 shape_type=5,
                 energy_range=35,
                 max_num=3,
                 split_components=True
                 ):
        self.shape_type = shape_type
        self.energy_range = energy_range
        self.max_num = max_num
        self.split_components = split_components

    def transform(self, iterable_seq=None):
        graphs = rnashapes_to_eden(iterable_seq, shape_type=self.shape_type, energy_range=self.energy_range,
                                   max_num=self.max_num, split_components=self.split_components)
        return graphs


class RNASynth():

    """Synthesizer class for RNA sequences. Produces new sequences similar in structure to the given sample set.

    Larger help explanation.
    Multi-line...

    Parameters
    ----------
    n_synthesized_sequences_per_seed_sequence : int (default 2)
            Option for setting the number of synthesized sequences per constraint.


    instance_score_threshold : int (default 0)
            Predicted score threshold for filtering synthesized sequences.


    shuffle_order : int (default 2)
            Eden.modifier.seq.seq_to_seq parameter.


    negative_shuffle_ratio : int (default 2)
            Number of negative sample sequences generated for each positive sample.


    vectorizer_complexity : int (default 2)
            eden.graph.Vectorizer parameter.


    max_num_graphs_per_seq: int (default 3)
            eden.converter.rna.rnashapes.rnashapes_to_eden parameter.
    """

    def __init__(self,
                 estimator=SGDClassifier(),
                 vectorizer=Vectorizer(),
                 pre_processor=PreProcessor(),
                 designer=AntaRNAv117Designer(),
                 constraint_extractor=ConstraintExtractor(),
                 n_synthesized_sequences_per_seed_sequence=3,
                 instance_score_threshold=0,
                 shuffle_order=2,
                 negative_shuffle_ratio=2,
                 n_jobs=-1,
                 cv=3,
                 n_iter_search=1
                 ):

        self.estimator = estimator
        self.vectorizer = vectorizer
        self.designer = designer
        self.pre_processor = pre_processor
        self.constraint_extractor = constraint_extractor

        self._n_synthesized_sequences_per_seed_sequence = n_synthesized_sequences_per_seed_sequence
        self._instance_score_threshold = instance_score_threshold
        self._shuffle_order = shuffle_order
        self._negative_shuffle_ratio = negative_shuffle_ratio
        self._n_jobs = n_jobs
        self._cv = cv
        self._n_iter_search = n_iter_search

        logger.debug('Instantiated an RNASynth object.')
        logger.debug(self.__dict__)

    def __repr__(self):
        obj_str = 'RNASynth:\n'
        obj_str += 'Dataset:\n'
        obj_str += 'shuffle_order: %d\n' % self._shuffle_order
        return obj_str

    def _binary_classification_setup(self, iterable_seq=None, negative_shuffle_ratio=None, shuffle_order=None):
        iter1, iter2 = tee(iterable_seq)
        iterable_graph = self.pre_processor.transform(iter1)
        iter3 = seq_to_seq(iter2, modifier=shuffle_modifier,
                           times=negative_shuffle_ratio, order=shuffle_order)
        iterable_graph_neg = self.pre_processor.transform(iter3)
        return iterable_graph, iterable_graph_neg

    def fit(self, iterable_seq):
        iterable_graph, iterable_graph_neg = self._binary_classification_setup(
            iterable_seq=iterable_seq,
            negative_shuffle_ratio=self._negative_shuffle_ratio,
            shuffle_order=self._shuffle_order)
        self.estimator = fit(iterable_graph,
                             iterable_graph_neg,
                             self.vectorizer,
                             n_jobs=self._n_jobs,
                             cv=self._cv,
                             n_iter_search=self._n_iter_search)
        return self

    def __design(self, iterable_graph):
        iterable_graph = self.vectorizer.annotate(
            iterable_graph, estimator=self.estimator)

        iterable = self.constraint_extractor.extract_constraints(
            iterable_graph)
        for (dot_bracket, seq_constraint, gc_content, fasta_id) in iterable:
            for count in range(self._n_synthesized_sequences_per_seed_sequence):
                sequence = self.designer.design(
                    (dot_bracket, seq_constraint, gc_content))
                header = fasta_id + '_' + str(count)
                yield header, sequence

    def __filter(self, iterable_seq):
        iter1, iter2 = tee(iterable_seq)
        iterable_graph = self.pre_processor.transform(iter1)
        predictions = self.vectorizer.predict(iterable_graph, self.estimator)

        for prediction, seq in izip(predictions, iter2):
            if prediction > self._instance_score_threshold:
                yield seq

    def sample(self, iterable_seq):
        iterable_graph = self.pre_processor.transform(iterable_seq)
        iterable_seq = self.__design(iterable_graph)
        iterable_seq = self.__filter(iterable_seq)

        return iterable_seq

    def fit_sample(self, iterable_seq):
        iterable_seq, iterable_seq_ = tee(iterable_seq)
        self.fit(iterable_seq)
        iterable_seq = self.sample(iterable_seq_)

        return iterable_seq


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to RNASynthesizer module.')

    rfam_id = 'RF01685'
    iterable_seq = fasta_to_sequence(
        'http://rfam.xfam.org/family/%s/alignment?acc=%s&format=fastau&download=0' % (rfam_id, rfam_id))
    initializer = Initializer()
    synthesizer = initializer.init_synthesizer()
    iter_seq = synthesizer.fit_sample(iterable_seq)
    for item in iter_seq:
        print item
