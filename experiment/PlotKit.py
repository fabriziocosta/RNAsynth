#!/usr/bin/env python

import os
from os.path import expanduser
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import pandas
from pandas import *
import argparse
import scipy as sp
from scipy.optimize import curve_fit
from numpy import exp


def get_column(dframe,col_name):
	
	listdf = dframe.describe().loc[col_name].tolist()
	
	return listdf


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


def cBox_plot(df1,df2,measure):
	"""
	Draws box plots of two sets of data in one figure.
	Accepts pandas dataframes as input.
	"""
	fig, ax1 = plt.subplots(figsize=(10,6))
	fig.canvas.set_window_title('BoxPlot Comparisons')
	plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
	
	plt.ylim(0.98,1.001)
	
	ax1.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
	
	ax1.set_title('%s BoxPlot Comparison' %measure)
	ax1.set_xlabel('Dataset Percentage Used for Training')
	ax1.set_ylabel('%s' %measure)
	
	bp1 = df1.boxplot(notch=0, sym='+', vert=1, whis=1.5)
	plt.setp(bp1['boxes'], color='red')
	plt.setp(bp1['whiskers'], color='red')
	plt.setp(bp1['medians'], linewidth=3, color='blue')
	plt.setp(bp1['fliers'], color='red', marker='+')
	
	bp2 = df2.boxplot(notch=0, sym='+', vert=1, whis=1.5)
	plt.setp(bp2['boxes'], color='green')
	plt.setp(bp2['whiskers'], color='green')
	plt.setp(bp2['medians'], linewidth=3, color='orange')
	plt.setp(bp2['fliers'], color='green', marker='+')
	
	
	plt.figtext(0.80, 0.06,  'True Samples' , \
		backgroundcolor='blue', color='white', weight='roman', \
		size='medium')
	plt.figtext(0.80, 0.025, 'Mixed Samples', \
		backgroundcolor='orange', \
		color='white', weight='roman', size='medium')
	plt.legend(loc = 4)
	plt.show()
	
	return True


if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--measure',required = True)
	parser.add_argument('-s', '--save', action = 'store_true')
	parser.add_argument('-f1', '-file1', required = True)
	parser.add_argument('-f2', '-file2', required = True)
	args = parser.parse_args()
	my_measure = args.measure
	my_image = my_measure + '.png'
	home = expanduser("~")
	my_path = os.path.join(home,my_image)
	
	file1 = args.file1
	file2 = args.file2
	
	x = np.array([10,20,30,40,50,60,70,80,90,100])

	dframe1 = read_csv(file1,header=0)
	dframe2 = read_csv(file2,header=0)
	
	yt = np.array(get_column(dframe1,'mean'))
	ys = np.array(get_column(dframe2,'mean'))
	
	#xpDecay_plot(x, yt, ys, my_measure)
	cBox_plot(dframe1, dframe2, my_measure)
	
	if args.save:
		plt.savefig(my_path)
		
	plt.show()

