#!/usr/bin/env python

import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import argparse
import scipy as sp
from scipy.optimize import curve_fit
from numpy import exp


def learning_curve_function(x, a, b):
	return a * (1 - np.exp(-b * x))

def draw_learning_curve(data_A = None , data_B = None, measure = None, x=None):
	"""
	Accepts as input an iterator over lists of numbers.
	Draws the exponential decay grpah over the means of lists.
	"""
	x=np.array(x)
	mean_originals = []
	for originals in data_A:
		mean_originals.append(mean(np.array(originals)))

	mean_originals_and_samples = []
	for originals_and_samples in data_B:
		mean_originals_and_samples.append(mean(np.array(originals_and_samples)))

	a, b = curve_fit(learning_curve_function, x, mean_originals)
	c, d = curve_fit(learning_curve_function, x, mean_originals_and_samples)

	x_fit = np.linspace(x.min(),x.max(),100)
	mean_originals_fit = learning_curve_function(x_fit, *a)
	mean_originals_and_samples_fit = learning_curve_function(x_fit, *c)

	fig, ax1 = plt.subplots(figsize=(10,6))
	fig.canvas.set_window_title('Exponential Decay Learning Curves')
	plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

	ax1.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
	ax1.set_title('Learning Curve Comparison for %s' %measure)
	ax1.set_xlabel('Dataset Percentage Used for Training')
	ax1.set_ylabel('%s Value' %measure)
	
	delta=0.25
	plt.boxplot(data_A,positions=x-delta)
	plt.plot(x, mean_originals, 'ro', label='')
	plt.plot(x_fit, mean_originals_fit, 'r-', label='True Samples')

	plt.boxplot(data_B,positions=x+delta)
	plt.plot(x, mean_originals_and_samples, 'go', label='')
	plt.plot(x_fit, mean_originals_and_samples_fit, 'g-', label='Mixed Samples')
	plt.grid()
	plt.legend(loc = 'lower right')
	plt.show()


if __name__ == "__main__":

	a_t = [[0.9369060577707512, 0.90133724392800463, 0.91931597267587639, 0.97606301975925636, 0.96563464634494833] , [0.9369060577707512, 0.90133724392800463, 0.91931597267587639, 0.97606301975925636, 0.96563464634494833],[0.9369060577707512, 0.90133724392800463, 0.91931597267587639, 0.97606301975925636, 0.96563464634494833]]
	a_s = [[0.96680658753051796, 0.94292599838634872, 0.95301928482157527, 0.9875969646313193, 0.98568403900205803],[0.96680658753051796, 0.94292599838634872, 0.95301928482157527, 0.9875969646313193, 0.98568403900205803],[0.96680658753051796, 0.94292599838634872, 0.95301928482157527, 0.9875969646313193, 0.98568403900205803]]

	draw_learning_curve(data_A = a_t , data_B = a_s , measure = 'ROC')