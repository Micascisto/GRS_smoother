# Pyhton script to smooth PDS Gamma Ray Spectrometer grids
# Original code by Michael Christoffersen, Brandon Tober
# Last edited 2019-11-22 by Stefano Nerozzi (stefano.nerozzi@utexas.edu)

import numpy as np
import sys
import math

# Constants
marsradius = 3396  # Mars radius in km
csize_in = 5 # size of input grid cells in degrees
csize_out = 0.5 # size of output grid cells in degrees
rad = 1.5  # search radius [sigma]
halfd = 265 # 50% signal radius [km]

avg_dist = 2

infile = 'cl_5x5.csv'
outfile = 'cl_5x5_gauss_1.5s.csv'
in_dir = '/home/stefano/Boulder_Halos/GRS/elements-v1/unsmoothed/5x5/'
out_dir = '/home/stefano/Boulder_Halos/GRS/elements-v1/processed/5x5/'
outfile = out_dir+outfile
infile = in_dir+infile

# calculate 1 sigma radius
sigma = halfd/0.67448975 # 50% is 0.67448975 sigma

# print some info
print('1-sigma radius: ' + str(int(sigma)) + ' km')
print('Search radius: ' + str(int(avg_dist*sigma)) + ' km')

# Haversine formula for distance
def havs(latp, lonp, lats, lons, marsradius):
    # latp, lonp are scalar
    # lats, lons are 2D arrays
    # Output is array with distance from (latp, lonp) to every point in (lats, lons)  
    lats = np.array(lats)
    lons = np.array(lons)
    dLat = np.radians(np.subtract(lats, latp))
    dLon = np.radians(np.subtract(lons, lonp))
    latp = np.radians(latp)
    lats = np.radians(lats)
 
    a = np.power(np.sin(dLat / 2),2) + np.cos(latp) * np.cos(lats) * np.power(np.sin(dLon / 2),2)
    c = 2 * np.arcsin(np.sqrt(a))
 
    return marsradius * c

def main():
  # Load input data file
  datalst = np.loadtxt(infile,delimiter=',')

  # Put input data into a grid
  latcells = 180/csize_in
  loncells = 360/csize_in
  datagrid = np.zeros((latcells,loncells))
  for k in range(len(datalst)):
    # Funky indexing to deal with order of points in file
    # This will break if ordering of points in input file changes
    # Probably better to actually read lat and lon and use that
    i = (latcells-1)-(k//loncells)
    j = k%loncells
    datagrid[i,j] = datalst[k,2]

  # Generate lat and lon arrays
  # lat 5x5 deg array looks like
  # [ 87.5, 87.5, 87.5, ... ]
  # [ 82.5, 82.5, 82.5, ... ]
  # [ 77.5, 77.5, 77.5, ... ]
  # [ ...., ...., ...., ... ]
  #
  # lon 5x5 deg array looks like
  # [ 2.5, 7.5, 12.5, ...]
  # [ 2.5, 7.5, 12.5, ...]
  # [ 2.5, 7.5, 12.5, ...]
  # [ ..., ..., ...., ...]
  #
  # These are used for haversine distance calculation
  lonl = np.matrix(np.linspace(csize_in/2, 360-csize_in/2, loncells))
  latl = np.matrix(np.linspace(90-csize_in/2, -90+csize_in/2, latcells))
  lon = np.dot(np.ones((datagrid.shape[0],1)), lonl)
  lat = np.dot(np.transpose(latl),np.ones((1,datagrid.shape[1])))

  # Make output grid
  outgrid = np.zeros((180/csize_out, 360/csize_out))
  olon = np.linspace(csize_out/2, 360-csize_out/2, 360/csize_out)
  olat = np.linspace(90-360/csize_out, -90+360/csize_out, 180/csize_out)

  # Do averaging, iterating over output grid
  for i in range(180/csize_out):
    for j in range(360/csize_out):
      # Haversine distance from each point to rest of grid
      dist = havs(olat[i], olon[j], lat, lon, marsradius)
      # Calculate weights for average based on gaussian function
      weights = math.e ** ((-dist**2) / (2*sigma**2))
      # Mask out points and weights outside search radius
      mask = dist < avg_dist*sigma
      subpts = datagrid[mask]
      print('weights: ' + str(weights))
      weights = weights[mask]
      print('weights masked: ' + str(weights))
      # Turn 9999.999 into NaN
      subpts[subpts == 9999.999] = np.NaN
      # Multiply values with their weights
      subpts = np.multiply(weights, subpts)
      # Calculate average ignoring NaNs
      outgrid[i,j] = np.nanmean(subpts)
      # Propagate errors (averaging improves errors by 1/sqrt(N averaged samples))
      #outgrid[i,j,1:] = outgrid[i,j,1:] / np.sqrt(np.nansum(weights))

  # Write out text file in similar format
  with open(outfile, "w") as f:
    for i in range(180/csize_out):
      for j in range(360/csize_out):
        f.write("{},{},{}\n".format(olat[i], olon[j], outgrid[i,j]))
   
main()