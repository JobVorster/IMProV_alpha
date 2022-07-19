import pandas as pd 
from glob import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
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
def monte_carlo_weighted_average(x,dx,weight,dweight,N):
	#x = 1/sum(weights)* sum(weight*x)
	weighted_average_arr = []
	for i in range(int(N)):
		rand_x = np.random.normal(x,dx)
		# rand_x = np.random.uniform(x-dx,x+dx)
		rand_weight = np.random.normal(weight,dweight)
		# rand_weight = np.random.uniform(weight-dweight,weight+dweight)
		weighted_average = sum(rand_x*rand_weight)/(sum(rand_weight))
		weighted_average_arr.append(weighted_average)
	return np.mean(weighted_average_arr),np.std(weighted_average_arr)*3
# data_dir = sys.argv[1]

posmaps = glob('./*.csv')
posmaps.sort()
feature_data = glob('./Features/*features.csv')
feature_data.sort()

for i,(posmap,feature) in enumerate(zip(posmaps,feature_data)):
	df_feature = pd.read_csv(feature,sep=',')
	df_posmap = pd.read_csv(posmap,sep=',')
	spot_RA = df_posmap['RA'].values
	spot_DEC = df_posmap['DEC'].values
	spot_RAerr = df_posmap['RAERR'].values
	spot_DECerr = df_posmap['DECERR'].values
	spot_FLUX = df_posmap['FLUX'].values
	spot_DFLUX = df_posmap['DFLUX'].values
	spot_VLSR = df_posmap['VLSR'].values

	feature_RA = []
	feature_RAERR = []
	feature_DEC = []
	feature_DECERR = []
	feature_VLSR = []

	for j in trange(len(df_feature.index)):
		xlim = [df_feature['RA min'].values[j],df_feature['RA max'].values[j]]
		ylim = [df_feature['Dec min'].values[j],df_feature['Dec max'].values[j]]
		inds = isolate(spot_RA,spot_DEC,xlim,ylim)
		subinds = np.where(spot_FLUX[inds]>0.5*df_feature['Peak'].values[j])
		RA_centroid,dRA_centroid = monte_carlo_weighted_average(
			spot_RA[inds][subinds],spot_RAerr[inds][subinds],spot_FLUX[inds][subinds],spot_DFLUX[inds][subinds],1e4)
		feature_RA.append(RA_centroid)
		feature_RAERR.append(dRA_centroid)
		DEC_centroid,dDEC_centroid = monte_carlo_weighted_average(
			spot_DEC[inds][subinds],spot_DECerr[inds][subinds],spot_FLUX[inds][subinds],spot_DFLUX[inds][subinds],1e4)
		feature_DEC.append(DEC_centroid)
		feature_DECERR.append(dDEC_centroid)
		feature_VLSR.append(df_feature['Center'].values[j])
	plt.scatter(feature_RA,feature_DEC,c = feature_VLSR,marker="$%d$"%(i+1),cmap='jet')
	df_processed = pd.DataFrame(columns=['RA','DEC','DRA','DDEC','VLSR','Peak','Peak err','Center','Center err','FWHM','FWHM err'])
	df_processed['RA'] = np.array(feature_RA)*1e-3
	df_processed['DEC'] = np.array(feature_DEC)*1e-3
	df_processed['DRA']= np.array(feature_RAERR)*1e-3
	df_processed['DDEC'] = np.array(feature_DECERR)*1e-3
	df_processed['VLSR'] = feature_VLSR
	for col in ['Peak','Peak err','Center','Center err','FWHM','FWHM err']:
		df_processed[col] = df_feature[col]
	df_processed.to_csv(posmap[:-4] + 'features.csv')
plt.show()	
