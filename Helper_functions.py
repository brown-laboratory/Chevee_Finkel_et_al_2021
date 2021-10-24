import numpy as np
import pandas as pd
import scipy as sp
from scipy import signal
from scipy import stats
from scipy.io import loadmat
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import math
from sklearn import metrics
from sklearn import linear_model
import pickle
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from matplotlib import cbook
from matplotlib import colors
from matplotlib.colors import Normalize


##############################################################################
def Resample_list(A):
    
    import numpy as np
    from random import sample

    Random_idx=[]
    for i in np.arange(len(A)): #this loops resamples the index of each unit within the list (with replacement)
        Random_idx.append(sample(np.arange(len(A)-1).tolist(),1)) # need to add -1 because of Zero indexing. the index of the last element of the list to sample from is len()-1
    Resampled_A=list( A[int(i)] for i in np.asarray(Random_idx))
    
    #the small recursive loop below is incorporated into Resample_list(A) to deal with the special case of only rows of one kind are selected when resampling for AUC.
    #In the specific case causing problems, the category is a 0 or a 1 at the end of each row, so to check that we don't only have 0s or 1s, we check the sums. If we do, we recall Resample_list(A)
    #This should not interfere with the normal function of Resample_list(A)
    X=np.asarray(Resampled_A)
    if (np.sum(X[:,-1]) == 0) | (np.sum(X[:,-1])==(np.shape(X)[0])):
       print(np.sum(X[:,-1]))
       Resampled_A=Resample_list(A)
    
    return Resampled_A

def  get_mean_and_95CI_resampled(A,iterations):
    
    import numpy as np
    from Helper_functions import Resample_list
    
    Samples_of_A=[]
    for iterations in np.arange(iterations): # loop trough the number of times you want
        Resampled_A=Resample_list(A) #call Resample_list to get random samples from A. Same size as A with replacement
        Samples_of_A.append(np.mean(Resampled_A, axis=0)) # append the mean of each iteration
    A_mean=np.mean(Samples_of_A, axis=0) # get the mean across all iterations
    A_sem=np.std(Samples_of_A, axis=0)#/np.sqrt(iterations) # get the SEM across all iterations
    A_95CI=1.96*A_sem
        #######NEED TO FIX THE CALCULATION OF SEM
    return A_mean, A_95CI

def Make_x_ysmooth(A, iterations, smooth=0.4):
    
    from scipy.interpolate import UnivariateSpline 
    from Helper_functions import get_mean_and_95CI_resampled
    import numpy as np
    A_mean, A_95CI = get_mean_and_95CI_resampled(A,iterations)
    x=np.arange(len(A_mean))
    spl_mean = UnivariateSpline(x, A_mean)
    spl_mean.set_smoothing_factor(smooth)
    y=spl_mean(x)
    
    return x,y, A_mean, A_95CI

    


class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
    
    
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
   