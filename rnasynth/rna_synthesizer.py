#!/usr/bin/env python


from itertools import tee
from itertools import izip
import logging

from sklearn.linear_model import SGDClassifier

from eden.converter.rna.rnafold import rnafold_to_eden
from eden.converter.fasta import fasta_to_sequence
from eden.converter.rna.rnashapes import rnashapes_to_eden
from eden.util import fit
from eden.graph import Vectorizer

from rnasynth.constraint_extractor import ConstraintExtractor
from rnasynth.rna_designer import AntaRNAv117Designer

from util.dataset import binary_classification_dataset_setup


logger = logging.getLogger(__name__)


def pre_process(iterable_seq, **opts):
    """
    DOCUMENTATION
    """
    graphs = rnashapes_to_eden(iterable_seq, shape_type=5, energy_range=35,
                               max_num=opts['max_num'], split_components=True)
    return graphs


class RNASynth(object):

    """Short one line.

    Larger help explanation.
    Multi-line...

    Parameters
    ----------
    n_synthesized_sequences_per_seed_sequence : int (default 2)
            Option for setting the number of synthesized sequences per constraint.


    instance_score_threshold : int (default 0)
            Predicted score threshold for filtering synthesized sequences.


    train_to_test_split_ratio : float (default 0.2)
            Ratio for splitting the sample dataset into train and test datasets.


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
                 designer=AntaRNAv117Designer(),
                 constraint_extractor=ConstraintExtractor(),
                 pre_processor=pre_process,
                 n_synthesized_sequences_per_seed_sequence=2,
                 instance_score_threshold=0,
                 train_to_test_split_ratio=0.2,
                 shuffle_order=2,
                 negative_shuffle_ratio=2,
                 vectorizer_complexity=2,
                 max_num_graphs_per_seq=3
                 ):

        self.estimator = estimator
        self.vectorizer = vectorizer
        self.designer = designer
        self.pre_processor = pre_processor
        self.constraint_extractor = constraint_extractor

        self._n_synthesized_sequences_per_seed_sequence = n_synthesized_sequences_per_seed_sequence
        self._instance_score_threshold = instance_score_threshold
        self._train_to_test_split_ratio = train_to_test_split_ratio
        self._shuffle_order = shuffle_order
        self._negative_shuffle_ratio = negative_shuffle_ratio
        self._vectorizer_complexity = vectorizer_complexity
        self._max_num_graphs_per_seq = max_num_graphs_per_seq

        logger.debug('Instantiated an RNASynth object.')
        logger.debug(self.__dict__)

    def __repr__(self):
        obj_str = 'RNASynth:\n'
        obj_str += 'Dataset:\n'
        obj_str += 'shuffle_order: %d\n' % self._shuffle_order
        return obj_str

    def fit(self, iterable_seq, n_iter_search=1):
        iterable_graph, iterable_graph_neg = binary_classification_dataset_setup(
            iterable_seq=iterable_seq,
            negative_shuffle_ratio=self._negative_shuffle_ratio,
            shuffle_order=self._shuffle_order)
        self.estimator = fit(iterable_graph,
                             iterable_graph_neg,
                             self.vectorizer,
                             n_jobs=-1,
                             cv=3,
                             n_iter_search=n_iter_search)
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
        iterable_graph = rnafold_to_eden(iter1)
        predictions = self.vectorizer.predict(iterable_graph, self.estimator)

        for prediction, seq in izip(predictions, iter2):
            if prediction > self._instance_score_threshold:
                yield seq

    def sample(self, iterable_seq):
        iterable_graph = rnafold_to_eden(iterable_seq)
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

    params = {'n_synthesized_sequences_per_seed_sequence': 2, 'instance_score_threshold': 0,
              'train_to_test_split_ratio': .2, 'shuffle_order': 2, 'negative_shuffle_ratio': 2,
              'vectorizer_complexity': 2, 'max_num_graphs_per_seq': 3}

    rfam_id = 'RF01685'
    iterable_seq = fasta_to_sequence(
        'http://rfam.xfam.org/family/%s/alignment?acc=%s&format=fastau&download=0' % (rfam_id, rfam_id))
    synthesizer = RNASynth(params)

    iter_seq = synthesizer.fit_sample(iterable_seq)
    for item in iter_seq:
        print item
