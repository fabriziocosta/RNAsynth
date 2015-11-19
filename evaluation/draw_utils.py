#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def learning_curve_function(x, a, b):
    return a * (1 - np.exp(-b * x))


def draw_learning_curve(data_a=None, data_b=None, measure=None, x=None, delta=0.1, scaling=100, fname=None):
    """
    Accepts as input an iterator over lists of numbers.
    Draws the exponential decay grpah over the means of lists.
    """
    x = np.array(x)
    mean_originals = []
    for originals in data_a:
        mean_originals.append(np.mean(np.array(originals)))

    mean_originals_and_samples = []
    for originals_and_samples in data_b:
        mean_originals_and_samples.append(
            np.mean(np.array(originals_and_samples)))

    a, b = curve_fit(learning_curve_function, x, mean_originals)
    c, d = curve_fit(learning_curve_function, x, mean_originals_and_samples)

    x_fit = np.linspace(x.min(), x.max(), 100)
    mean_originals_fit = learning_curve_function(x_fit, *a)
    mean_originals_and_samples_fit = learning_curve_function(x_fit, *c)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.canvas.set_window_title('Exponential Decay Learning Curves')
    # plt.subplots_adjust(left=0.04, right=0.35, top=0.9, bottom=0.25)

    ax1.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax1.set_title('Learning Curve Comparison for %s' % measure)
    ax1.set_xlabel('Dataset Percentage Used for Training')
    ax1.set_ylabel('%s Value' % measure)

    plt.boxplot(data_a, positions=(x + delta) * scaling, notch=False)
    plt.plot((x + delta) * scaling, mean_originals, 'ro', label='')
    plt.plot((x_fit) * scaling, mean_originals_fit, 'r-', label='Original')

    plt.boxplot(data_b, positions=(x - delta) * scaling, notch=False)
    plt.plot((x - delta) * scaling, mean_originals_and_samples, 'go', label='')
    plt.plot((x_fit) * scaling, mean_originals_and_samples_fit,
             'g-', label='Original+sampled')
    plt.grid()
    plt.legend(loc='lower right')
    if fname is not None:
        plt.savefig(fname)
    else:
        plt.show()


if __name__ == "__main__":

    prefix = 'RF00005'
    delta = 0.005
    scaling = 50

    a_t = [
        [0.61882589177433867, 0.86558548786434841, 0.81731460783511856, 0.87284584856326064, 0.69121999649044663],
        [0.9069060577707512, 0.91133724392800463, 0.92931597267587639,
            0.93606301975925636, 0.95563464634494833],
        [0.89369060577707512, 0.90133724392800463, 0.91931597267587639,
            0.95606301975925636, 0.96563464634494833]]
    a_s = [[0.63349937846848603, 0.92822770048379799, 0.74499106001252102, 0.91760562533385914, 0.73954523834430574],
           [0.96680658753051796, 0.94292599838634872, 0.95301928482157527,
               0.9875969646313193, 0.98568403900205803],
           [0.96680658753051796, 0.94292599838634872, 0.95301928482157527, 0.9875969646313193, 0.98568403900205803]]

    data_fractions = list(np.linspace(0.2, 1.0, 3))
    draw_learning_curve(data_a=a_t, data_b=a_s, x=data_fractions, measure='ROC',
                        delta=delta, scaling=scaling, fname=prefix + 'roc.pdf')
