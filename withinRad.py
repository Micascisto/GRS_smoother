import math
import numpy as np
import sys
import time

"""
extract points from csv that are within radius from desired point
"""

### INPUTS ###
infile = 'cl_5x5.csv'
outfile = 'cl_5x5_gauss_2s.csv'
in_dir = '/home/stefano/Boulder_Halos/GRS/elements-v1/unsmoothed/5x5/'
out_dir = '/home/stefano/Boulder_Halos/GRS/elements-v1/processed/5x5/'
outfile = out_dir+outfile
infile = in_dir+infile

rad = 2  # search radius [sigma]
halfd = 265 # 50% signal radius [km]
step = 0.5 # [degrees]
nd = 9999.999 # null value

### CODE ###
t0 = time.time()

# calculate the haversine distance from np array of points
def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a)) 
    r = 3389.5 # Radius of body in km (mean Mars radius)
    return c * r

# read csv to np array
data = np.genfromtxt(infile, delimiter=',')

# replace nd value with NaN
data[data==nd] = np.nan

# create masked array for NaN
masked_data = np.ma.masked_array(data, np.isnan(data))

# empty distance array to hold haversine distance between points in array
dist = np.zeros(len(data))

# empty weights array
weight = np.zeros(len(data))

# get len of out data based on step size
num_its = len(np.arange(-90,90 + step, step)) * len(np.arange(0,360 + step, step))

# output data array
out = np.zeros((num_its,masked_data.shape[1]))

# calculate 1 sigma radius
sigma = halfd/0.67448975 # 50% is 0.67448975 sigma

# print some info
print('1-sigma radius: ' + str(round(sigma,0)) + ' km')
print('Search radius: ' + str(round(rad*sigma,0)) + ' km')

# iterate through lat and lon range
count = 0
for i in np.arange(-90,90 + step, step):
    for j in np.arange(0,360 + step, step):
        for k in range(len(dist)):
            dist[k] = haversine(i, j, data[k,0], data[k,1])
            weight[k] = math.e ** ((-dist[k]**2) / (2*sigma**2))
        # for this lat and lon, find distances within search radius
        idx = np.where(dist <= rad*sigma)
        # assign lat and lon to out array
        out[count,0] = i
        out[count,1] = j
        # calculate weighted average
        out[count,2:] = np.ma.filled(np.ma.average(masked_data[idx,2:], axis=1, weights=weight[idx]), np.nan)
        # propagate errors
        out[count,3] = out[count,3] / np.sqrt(np.nansum(weight[idx]))
        out[count,4] = out[count,4] / np.sqrt(np.nansum(weight[idx]))
        # increase iterator
        count += 1

# remove nan points
out = out[~np.isnan(out[:,2]),:]

# save csv
np.savetxt(outfile, out, delimiter=',', newline = '\n', comments = '', fmt='%.8f')

t1 = time.time()
print('Total Runtime: ' + str(round((t1 - t0),4)) + ' seconds')
