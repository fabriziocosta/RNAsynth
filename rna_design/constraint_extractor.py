#!/usr/bin/env python

import re
import logging
import networkx as nx


def dict_to_string(dictionary):
	"""
	Generic function to build a sequential string of dictionary values.
	"""
	st = ''
	for i in range(len(dictionary)):
		st = st + dictionary[i]
		
	return st


def build_nodes_dict(g):
	"""
	Builds a dictionary of key = node value = nucleotide out of the graph.
	"""
	nodes_dict = {}
	for node, data in g.nodes_iter(data=True):
		nodes_dict.update({node:data['label']})
		
	return nodes_dict


def build_generic_nodes_dict(g, padding = 'A'):
	"""
	Builds a dictionary of key = node value = padding out of the graph.
	"""
	nodes_dict = {}
	for node, data in g.nodes_iter(data=True):
		nodes_dict.update({node:padding})
		
	return nodes_dict


def compute_gc_content(g):
	"""
	Function to calculate the GC content of all subgraphs in a graph set.
	"""
	gc_content = 0
	for node, data in g.nodes_iter(data=True):
		if (data['label']== 'G') or (data['label']== 'C'):
			gc_content += 1
	gc_content = float(gc_content) / float(nx.number_of_nodes(g))
	
	return gc_content


def get_basepair_list(g):
	"""
	Accepts single graph as input.
	Returns a list of all base pairs in the folding.
	"""
	list_bpairs = []
	for line in nx.generate_edgelist(g):
		if line.find('basepair') > 0:
			list_bpairs.append((int(line.split(' ')[0]),int(line.split(' ')[1])))
			
	return list_bpairs


def importance_based_graph_cut(g, threshold):
	"""
	Removes nodes with importance below the threshold from g.
	"""
	for node, data in g.nodes_iter(data=True):
		if float(data['importance']) < threshold:
			g.remove_node(node)
			
	return


def get_importance_list(g, threshold, adjacency, importance = 1):
	"""
	Generates a list of important nodes in a graph.
	Importance is based on the importance number being greater than threshold,
	and adjaceny factor being greater than or equal to radius.
	Returns the complement list if importance=-1.
	"""
	graph = g.copy()
	nodes_list = []
	importance_based_graph_cut(graph, threshold)
	for component in nx.connected_components(graph):
		if len(component) >= adjacency:
			nodes_list = nodes_list + component
	if importance == 1:
		importance_list = nodes_list
	elif importance == -1:
		importance_list = [node for node in g.nodes() if node not in nodes_list]
		
	return importance_list


def find_paired_nodes(g):
	"""
	Returns a list containing all paired nodes in a graph.
	"""
	paired_list = []
	for line in nx.generate_edgelist(g):
		if ('basepair' in line):
			if not (int(line.split(' ')[0]) in paired_list):
				paired_list.append(int(line.split(' ')[0]))
			if not (int(line.split(' ')[1]) in paired_list):
				paired_list.append(int(line.split(' ')[1]))
	
	return paired_list


def pair_based_graph_cut(g):
	"""
	Removes paired nodes from g.
	"""
	for node in find_paired_nodes(g):
		g.remove_node(node)
	
	return


def find_unpaired_regions(g, adjacency):
	"""
	Generates a list of unpaired nodes in a graph.
	and adjaceny factor being greater than or equal to radius.
	Returns a list containing regions in which the number of unpaired nodes
	is greater than "adjacency".
	"""
	graph = g.copy()
	unpaired_nodes_list = []
	pair_based_graph_cut(graph)
	for component in nx.connected_components(graph):
		if len(component) >= adjacency:
			unpaired_nodes_list = unpaired_nodes_list + component

	return unpaired_nodes_list


def extract_constraint_string(g, threshold, adjacency = 1 , padding= 'N'):
	"""
	Generates the anatRNA input constraint string based on Networkx graph representation.
	Adjacent nodes above the threshold show up in the output string.
	"""
	Cstr_dict = build_nodes_dict(g)
	for node in get_importance_list(g, threshold, adjacency,importance=-1):
		Cstr_dict[node] = padding

	return dict_to_string(Cstr_dict)


def extract_dot_struct(g,threshold, importance_adjacency = 1, unpaired_adjacency = 1):
	"""
	Generates the anatRNA input dot-bracket notation constraint based on networkx graph representation.
	Base pairs above the importance threshold appear in the output string.
	"""
	list_bpairs = get_basepair_list(g)
	list_unpaired = find_unpaired_regions(g, unpaired_adjacency)
	dic_dot_str = build_generic_nodes_dict(g)
	importance_list = get_importance_list(g, threshold, importance_adjacency)
	dot_list = []
	
	for i,j in list_bpairs:
		if i in importance_list and j in importance_list:
			dic_dot_str[i] = '('
			dic_dot_str[j] = ')'
	for unpaired_node in list_unpaired:
		dic_dot_str[unpaired_node] = '.'

	return dict_to_string(dic_dot_str)


def extract_constraints(graphs, Cseq_threshold, Cseq_adjacency, dotnot_threshold, dotnot_adjacency, unpaired_adjacency):
	"""
	Generator function responsible for composing anatRNA dot-bracket and constraint string.
	Accepts different threshold and adjacency values for composing anatRNA dot-bracket and constraint string.
	"""
	for g in graphs:
		fasta_id = g.graph['id']
		gc_content = compute_gc_content(g)
		Cseq = extract_constraint_string(g, Cseq_threshold, Cseq_adjacency)
		struct = extract_dot_struct(g,dotnot_threshold, dotnot_adjacency, unpaired_adjacency)
		
		yield struct,Cseq,gc_content,fasta_id


if __name__ == "__main__":

	logging.basicConfig(level=logging.INFO)
	logger = logging.getLogger(__name__)
	logger.info('Call to constraint_extractor package.')
