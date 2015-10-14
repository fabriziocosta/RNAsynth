#!/usr/bin/env python

import logging
import networkx as nx


class ConstraintExtractor():

    def __init__(self,
                 importance_threshold_sequence_constraint=0,
                 min_size_connected_component_sequence_constraint=1,
                 importance_threshold_structure_constraint=0,
                 min_size_connected_component_structure_constraint=1,
                 min_size_connected_component_unpaired_structure_constraint=1
                 ):
        self.importance_threshold_sequence_constraint = importance_threshold_sequence_constraint
        self.min_size_connected_component_sequence_constraint = min_size_connected_component_sequence_constraint
        self.importance_threshold_structure_constraint = importance_threshold_structure_constraint
        self.min_size_connected_component_structure_constraint = min_size_connected_component_structure_constraint
        self.min_size_connected_component_unpaired_structure_constraint = min_size_connected_component_unpaired_structure_constraint

    def extract_constraints(self, graphs):
        """
        Generator function which yields sequence and structure constraint strings extracted from an annotated Networkx graph.
        Accepts connectivity values and thresholds for sequence and structure constraints for the same graph separately .
        """
        for g in graphs:
            fasta_id = g.graph['id']
            gc_content = self._compute_gc_content(g)
            cseq = self._extract_constraint_string(
                g, self.importance_threshold_sequence_constraint, self.min_size_connected_component_sequence_constraint)
            struct = self._extract_dot_struct(
                g, self.importance_threshold_structure_constraint, self.min_size_connected_component_structure_constraint,
                self.min_size_connected_component_unpaired_structure_constraint)

            yield struct, cseq, gc_content, fasta_id

    def _extract_constraint_string(self,
                                   graph,
                                   importance_threshold_sequence_constraint,
                                   min_size_connected_component_sequence_constraint,
                                   padding='N'):
        """
        Generates a sequence constraint string from an annotated Networkx graph.
        Adjacent nodes with the connectivity above the threshold show up
        in the output string as actual nucleotides.
        Other nodes appear as padding in the output string.
        """
        cstr_dict = self._build_nodes_dict(graph)
        for node in self._get_importance_list(graph, importance_threshold_sequence_constraint,
                                              min_size_connected_component_sequence_constraint, importance=-1):
            cstr_dict[node] = padding
        cstr = self._dict_to_string(cstr_dict)
        return cstr

    def _extract_dot_struct(self,
                            graph,
                            importance_threshold_structure_constraint,
                            min_size_connected_component_structure_constraint,
                            min_size_connected_component_unpaired_structure_constraint):
        """
        Generates a dot-bracket structure constraint string from an annotated Networkx graph.
        Base pairs above the importance threshold appear in the output string.
        """
        list_bpairs = self._get_basepair_list(graph)
        list_unpaired = self._find_unpaired_regions(
            graph, min_size_connected_component_unpaired_structure_constraint)
        dic_dot_str = self._build_generic_nodes_dict(graph)
        importance_list = self._get_importance_list(
            graph, importance_threshold_structure_constraint, min_size_connected_component_structure_constraint)

        for i, j in list_bpairs:
            if i in importance_list and j in importance_list:
                dic_dot_str[i] = '('
                dic_dot_str[j] = ')'
        for unpaired_node in list_unpaired:
            dic_dot_str[unpaired_node] = '.'
        cstruct = self._dict_to_string(dic_dot_str)
        return cstruct

    def _dict_to_string(self, dictionary):
        """
        Generic function to build a sequential string of dictionary values.
        """
        st = ''
        for i in range(len(dictionary)):
            st = st + dictionary[i]
        return st

    def _build_nodes_dict(self, graph):
        """
        Builds a dictionary of key = node value = nucleotide out of the graph.
        """
        nodes_dict = {}
        for node, data in graph.nodes_iter(data=True):
            nodes_dict.update({node: data['label']})
        return nodes_dict

    def _build_generic_nodes_dict(self, graph, padding='A'):
        """
        Builds a dictionary of key = node value = padding out of the graph.
        """
        nodes_dict = {}
        for node, data in graph.nodes_iter(data=True):
            nodes_dict.update({node: padding})
        return nodes_dict

    def _compute_gc_content(self, graph):
        """
        Function to calculate the GC content of all subgraphs in a graph set.
        """
        gc_content = 0
        for node, data in graph.nodes_iter(data=True):
            if (data['label'] == 'G') or (data['label'] == 'C'):
                gc_content += 1
        gc_content = float(gc_content) / float(nx.number_of_nodes(graph))
        return gc_content

    def _get_basepair_list(self, graph):
        """
        Accepts single graph as input.
        Returns a list of all base pairs in the folding.
        """
        list_bpairs = []
        for line in nx.generate_edgelist(graph):
            if line.find('basepair') > 0:
                list_bpairs.append(
                    (int(line.split(' ')[0]), int(line.split(' ')[1])))
        return list_bpairs

    def _importance_based_graph_cut(self, graph, threshold):
        """
        Removes nodes with importance below the threshold from g.
        """
        for node, data in graph.nodes_iter(data=True):
            if float(data['importance']) < threshold:
                graph.remove_node(node)
        return

    def _get_importance_list(self, graph, threshold, adjacency, importance=1):
        """
        Generates a list of important nodes in a graph.
        Importance is based on the importance number being greater than threshold,
        and adjaceny factor being greater than or equal to radius.
        Returns the complement list if importance=-1.
        """
        graph_c = graph.copy()
        nodes_list = []
        self._importance_based_graph_cut(graph_c, threshold)
        for component in nx.connected_components(graph_c):
            if len(component) >= adjacency:
                nodes_list = nodes_list + component
        if importance == 1:
            importance_list = nodes_list
        elif importance == -1:
            importance_list = [
                node for node in graph.nodes() if node not in nodes_list]
        return importance_list

    def _find_paired_nodes(self, graph):
        """
        Returns a list containing all paired nodes in a graph.
        """
        paired_list = []
        for line in nx.generate_edgelist(graph):
            if ('basepair' in line):
                if not (int(line.split(' ')[0]) in paired_list):
                    paired_list.append(int(line.split(' ')[0]))
                if not (int(line.split(' ')[1]) in paired_list):
                    paired_list.append(int(line.split(' ')[1]))
        return paired_list

    def _pair_based_graph_cut(self, graph):
        """
        Removes paired nodes from graph.
        """
        for node in self._find_paired_nodes(graph):
            graph.remove_node(node)
        return

    def _find_unpaired_regions(self, graph, adjacency):
        """
        Generates a list of unpaired nodes in a graph.
        and adjaceny factor being greater than or equal to radius.
        Returns a list containing regions in which the number of unpaired nodes
        is greater than "adjacency".
        """
        graph_c = graph.copy()
        unpaired_nodes_list = []
        self._pair_based_graph_cut(graph_c)
        for component in nx.connected_components(graph_c):
            if len(component) >= adjacency:
                unpaired_nodes_list = unpaired_nodes_list + component
        return unpaired_nodes_list


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    logger.info('Call to constraint_extractor package.')
    CE = ConstraintExtractor()
