#!/usr/bin/env python

import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import argparse
import scipy as sp
from scipy.optimize import curve_fit
from numpy import exp


def func(x, a, b):
	return a * (1 - np.exp(-b * x))


def xpDecay_plot(x, yt, ys, measure):
	"""
	Plots the exponential decay curve of two sets of datasets in one figure.
	Accepts numpy arrays as input.
	"""
	a, b = curve_fit(func, x, yt)
	c, d = curve_fit(func, x, ys)

	xnew = np.linspace(x.min(),x.max(),100)
	ytnew = func(xnew, *a)
	ysnew = func(xnew, *c)

	fig, ax1 = plt.subplots(figsize=(10,6))
	fig.canvas.set_window_title('Exponential Decay Learning Curves')
	plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

	ax1.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
	ax1.set_title('Learning Curve Comparison for %s' %measure)
	ax1.set_xlabel('Dataset Percentage Used for Training')
	ax1.set_ylabel('%s Value' %measure)
	
	plt.plot(xnew, ytnew, 'r-', label='True Samples')
	plt.plot(xnew, ysnew, 'g-', label='Mixed Samples')
	plt.legend(loc = 4)
	plt.show()
	return True


def xpDecay(iterable_t = None , iterable_s = None, measure = None):
	"""
	Accepts as input an iterator over lists of numbers.
	Draws the exponential decay grpah over the means of lists.
	"""
	means_t = []
	for list in iterable_t:
		means_t.append(mean(np.array(list)))
	print means_t

	means_s = []
	for list in iterable_s:
		means_s.append(mean(np.array(list)))
	print means_s
	x = np.array([n for n in range(10, (len(means_t) + 1)*10, 10)])
	print x
	xpDecay_plot(x , means_t , means_s , measure)


if __name__ == "__main__":

	a_t = [[0.9369060577707512, 0.90133724392800463, 0.91931597267587639, 0.97606301975925636, 0.96563464634494833] , [0.9369060577707512, 0.90133724392800463, 0.91931597267587639, 0.97606301975925636, 0.96563464634494833],[0.9369060577707512, 0.90133724392800463, 0.91931597267587639, 0.97606301975925636, 0.96563464634494833]]
	a_s = [[0.96680658753051796, 0.94292599838634872, 0.95301928482157527, 0.9875969646313193, 0.98568403900205803],[0.96680658753051796, 0.94292599838634872, 0.95301928482157527, 0.9875969646313193, 0.98568403900205803],[0.96680658753051796, 0.94292599838634872, 0.95301928482157527, 0.9875969646313193, 0.98568403900205803]]

	xpDecay(iterable_t = a_t , iterable_s = a_s , measure = 'TEST')

