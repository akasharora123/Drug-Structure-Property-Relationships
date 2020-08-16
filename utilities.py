#!/usr/local/bin/env python
# -*- coding: utf-8 -*-

"""
A few utility functions:

-- polymerize: polymerize a molecule, i.e.,
   translate SMILES by adding repeat units in between 

-- inverse_transform_labels: To tranform the reposnse (one-hot-encoded) 
   variables back into label names for a multi-class multi-label problem

-- setup_plot: functiont to customize plot

"""

import os
import math
import numpy as np
import scipy as sp
import pandas as pd
import inspect
import rdkit
from rdkit import Chem
from rdkit.Chem import RDKFingerprint
import rdkit.Chem.Descriptors as Descriptors
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import PandasTools
from sklearn import metrics
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


def polymerize(Smiles, repeats):

    start = Smiles[0]+'6'+Smiles[1:]
    end = Smiles+'6'
    return start+repeats*Smiles+end


def setup_plot(ax, lgflag):
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(axis='both', which='major', direction='in',
                    top=True, right=True, length=7, width=0.25)
    ax.tick_params(axis='both', which='minor', direction='in',
                  top=True, right=True, length=3.5, width=0.25)
    if(lgflag=='Y'):
       lg = ax.legend()
       lg.get_frame().set_linewidth(0.25)


def compute_metrics(y_true, y_pred):
    
    mse = metrics.mean_squared_error(y_true,y_pred)
    r2 = metrics.r2_score(y_true,y_pred)    
    return math.sqrt(mse), r2


def plot_reg(ax, title, prop, y_train, y_train_pred, y_test=[], y_test_pred=[]):

    """Scatter plot of the predicted vs true targets."""
    ax.plot([y_train.min(), y_train.max()],
            [y_train.min(), y_train.max()],
            '--r', linewidth=2)
   
    if(any(y_test)): 
        ax.scatter(y_train, y_train_pred, alpha=0.2, label='train')
        rmse, r2 = compute_metrics(y_train, y_train_pred)
        ax.text(0.05, 0.88, 'Train RMSE: {:0.2f}, Train $R^2$: {:0.2f}'.format(rmse, r2), 
                            transform=ax.transAxes, fontsize=12)
        ax.scatter(y_test, y_test_pred, alpha=0.2, label='test')
        rmse, r2 = compute_metrics(y_test, y_test_pred)
        ax.text(0.05, 0.80, 'Test RMSE: {:0.2f}, Test $R^2$: {:0.2f}'.format(rmse, r2), 
                            transform=ax.transAxes, fontsize=12)
        
    else:
        ax.scatter(y_train, y_train_pred, alpha=0.2, label='test+train')
        rmse, r2 = compute_metrics(y_train, y_train_pred)
        ax.text(0.05, 0.88, 'RMSE$=${:0.2f}, $R^2=${:0.2f}'.format(rmse, r2), 
                            transform=ax.transAxes, fontsize=12)
        
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(axis='both', which='major', direction='in', 
                    top=True, right=True, length=7, width=0.25)
    ax.tick_params(axis='both', which='minor', direction='in', 
                  top=True, right=True, length=3.5, width=0.25)
    lg = ax.legend(loc=4)
    lg.get_frame().set_linewidth(0.25)
    ax.set_title(title, fontsize=16)
    ax.set_xlabel('Actual '+prop, fontsize=16)
    ax.set_ylabel('Predicted '+prop, fontsize=16) 



def plot_style():

    dark_gray = ".15"
    light_gray = ".8"
        
    style_dict = { 
                "figure.facecolor": "white",
                "text.color": dark_gray,
                "axes.labelcolor": dark_gray,
                "legend.frameon": True,
                "legend.numpoints": 1,
                "legend.scatterpoints": 1,
                "legend.fancybox": False,
                "legend.labelspacing": 0.40,
                "legend.handlelength": 1.25,
                "legend.handletextpad": 0.40,
                "legend.borderaxespad": 0.75,
                "legend.borderpad": 0.40,
                "xtick.direction": "out",
                "ytick.direction": "out",
                "xtick.color": dark_gray,
                "ytick.color": dark_gray,
                "axes.axisbelow": True,
                "image.cmap": "Greys",
                "font.family": ["sans-serif"],
                "font.sans-serif": [
                    "Arial",
                    "Liberation Sans",
                    "Bitstream Vera Sans",
                    "sans-serif",
                ],  
                "grid.linestyle": "-",
                "axes.grid": True,
                "lines.solid_capstyle": "round",
                "axes.facecolor": "white",
                "axes.edgecolor": light_gray,
                "axes.linewidth": 1.25,
                "grid.color": light_gray,
                "xtick.major.size": 0,
                "ytick.major.size": 0,
                "xtick.minor.size": 0,
                "ytick.minor.size": 0,
                "text.usetex": True,
                "text.latex.preamble": r'\usepackage{amsmath}'
            }

    return style_dict
