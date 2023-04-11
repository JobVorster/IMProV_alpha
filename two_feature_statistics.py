"""
IMPRoV 0.5.0

This is the second program of Interactive Proper motions for VLBI (IMProV).

AIM OF THE PROGRAM
The aim of this program is to calculate feature statistics of detected maser spots on VLBI images. The program requires the user to input the foldernames containing the regions of interest (ROIs) as well as the foldername containing the original data. The statistics of the detected maser spots within the ROIs are then calculated, including their flux weighted centroid and a Gaussian fit is done to the data. If the Gaussian fit does not succeed, the parameters are replaced by 999.

Author: Job Vorster
Date: April 11, 2023

Usage: python two_feature_statistics.py SPOTMAPS_FOLDERNAME FEATURE_BOXES_FOLDERNAME

Requirements:

    Python 3
    numpy
    pandas
    matplotlib
    scipy
    skspatial
    glob

Any questions, comments, or suggestions can be sent to:
jobvorster8@gmail.com
"""

#########################################################################################
#                                                                                       #
#                               Importing Libraries                                     #
#                                                                                       #
#########################################################################################


from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
from skspatial.objects import Line
import sys
from scipy.optimize import curve_fit


#########################################################################################
#                                                                                       #
#                               Function Definitions                                    #
#                                                                                       #
#########################################################################################



def isolate(x,y,xlim,ylim):
    '''Gives indices for elements in x,y range.
    Parameters :
    ------------
    x : array-like
    x coordinates of element set.
    y : array-like
    y coordinates of element set.
    xlim : array-like of form :[xlimmin,xlimmax]
    specify maximum and minimum x values for elements to isolate.
    ylim : array-like of form :[ylimmin,ylimmax]
    specify maximum and minimum y values for elements to isolate.
    
    Returns :
    --------
    ind : array-like
    Array containing indices of isolated elements. Can then be easily called, e.g. x-coords : x[ind]'''
    if len(xlim) == 0 or len(ylim) ==0:
        return []
    else:
        ind = np.where(np.logical_and(np.logical_and(y>= ylim[0],y<=ylim[1]),np.logical_and(x <= xlim[1],x>=xlim[0])))
        return ind[0]

def lin_func(x,slope,intercept):
	return x*slope + intercept

def gaussian(x,A,mu,sigma):
	return A*np.exp(-(x-mu)**2/(2*sigma**2))

def flux_weighted_centroid(RA,DEC,FLUX):
	return [(np.array(x)*np.array(FLUX))/sum(FLUX) for x in [RA, DEC]]

def run_mc(RA,DEC,RAerr,DECerr,FLUX,FLUXerr,N):
	RA = RA.values
	DEC = DEC.values
	RAerr = RAerr.values
	DECerr = DECerr.values
	FLUX = FLUX.values
	FLUXerr = FLUXerr.values	
	inds = np.where(FLUX > 0.5*max(FLUX))[0]
	print('Inds = %s, if it is none. There is a problem'%(str(inds)))


	RA_mc = [np.random.normal(ra,raerr,size=N) for (ra,raerr) in zip(RA[inds],RAerr[inds])]
	DEC_mc = [np.random.normal(dec,decerr,size=N) for (dec,decerr) in zip(DEC[inds],DECerr[inds])]
	FLUX_mc = [np.random.normal(f,ferr,size=N) for (f,ferr) in zip(FLUX[inds],FLUXerr[inds])]

	ra_c = np.zeros(N)
	dec_c = np.zeros(N)

	for i in range(len(inds)):
		ra_c += (np.array(RA_mc)[i,:]*np.array(FLUX_mc)[i,:])/np.sum(FLUX_mc,axis=0)
		dec_c += (np.array(DEC_mc)[i,:]*np.array(FLUX_mc)[i,:])/np.sum(FLUX_mc,axis=0)

	print("The shape of ra_c is %s. If it is equal to %d it is correct."%(str(np.shape(ra_c)),N))
	centre_ra = np.mean(ra_c)
	std_ra = 3*np.std(ra_c)
	centre_dec = np.mean(dec_c)
	std_dec = np.std(dec_c)
	return centre_ra,std_ra,centre_dec,std_dec


#########################################################################################
#                                                                                       #
#                               Main code 			                                    #
#                                                                                       #
#########################################################################################


if __name__ == "__main__":

	#The folder should contain only the spotmaps or boxes. They should be named with EpXXSOURCE_NAME.csv with XX the epoch number and SOURCE_NAME the name of the astrophysical source.
	spotmaps_foldername = sys.argv[1]
	boxes_foldername = sys.argv[2]


	spotmaps_f = glob('%s/*.csv'%(spotmaps_foldername))
	spotmaps_f.sort()
	boxes_f = glob('%s/*.csv'%(boxes_foldername))
	boxes_f.sort()
	
	mc_iterations = 10000
	for i,(spots,box) in enumerate(zip(spotmaps_f,boxes_f)):
		feature_pars = {}
		df_spots = pd.read_csv(spots)
		df_box = pd.read_csv(box)

		for init_pars in ['RA','DRA','DEC','DDEC','Peak','Peak err','Center','Center err','FWHM','FWHM err','VLSR']:
			feature_pars[init_pars] = []

		for feature in range(len(df_box["RA min"])):

			xlim = [df_box['RA min'][feature],df_box['RA max'][feature]]
			ylim = [df_box['Dec min'][feature],df_box['Dec max'][feature]]
			inds = isolate(df_spots['RA'],df_spots['DEC'],xlim,ylim)

			try:
				inds_flux = df_spots['FLUX'].values[inds] > 0.5*max(df_spots['FLUX'].values[inds])
				p0 = [max(df_spots['FLUX'].values[inds][inds_flux]),df_spots['VLSR'].values[inds][inds_flux][np.argmax(df_spots['FLUX'].values[inds][inds_flux])],0.5]
				popt,pcov = curve_fit(gaussian,df_spots['VLSR'].values[inds][inds_flux],df_spots['FLUX'].values[inds][inds_flux],sigma = df_spots['DFLUX'].values[inds][inds_flux],absolute_sigma=True,p0=p0)
				feature_pars['Peak'].append(np.round(popt[0]))
				feature_pars['Center'].append(np.round(popt[1]))
				feature_pars['VLSR'].append(np.round(popt[1]))
				feature_pars['FWHM'].append(np.round(popt[2]*2*np.sqrt(np.log(2))))

				feature_pars['Peak err'].append(np.round(np.sqrt(np.diag(pcov))[0]))
				feature_pars['Center err'].append(np.round(np.sqrt(np.diag(pcov))[1]))
				feature_pars['FWHM err'].append(np.round(np.sqrt(np.diag(pcov))[2]))
			except:
				feature_pars['Peak'].append(999)
				feature_pars['Center'].append(999)
				feature_pars['VLSR'].append(999)
				feature_pars['FWHM'].append(999)

				feature_pars['Peak err'].append(999)
				feature_pars['Center err'].append(999)
				feature_pars['FWHM err'].append(999)

			print('Calculating Parameters for Feature %d of %d for Epoch %d of %d. N = %d'%(feature+1,len(df_box['RA min']),i+1,len(spotmaps_f),mc_iterations))
			mc_pars_mean = ['RA','DEC']
			mc_pars_std = ['DRA','DDEC']

			centre_ra,std_ra,centre_dec,std_dec = run_mc(df_spots['RA'][inds],df_spots['DEC'][inds],df_spots['RAERR'][inds],df_spots['DECERR'][inds],df_spots['FLUX'][inds],df_spots['DFLUX'][inds],mc_iterations)
			feature_pars['RA'].append(np.round(centre_ra))
			feature_pars['DRA'].append(np.round(std_ra))
			feature_pars['DEC'].append(np.round(centre_dec))
			feature_pars['DDEC'].append(np.round(std_dec))
		df_feature_pars = pd.DataFrame.from_dict(feature_pars)
		df_feature_pars.to_csv('Ep%dfeature_parameters.csv'%(i+1))