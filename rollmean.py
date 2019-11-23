# Pyhton script to smooth PDS Gamma Ray Spectrometer grids
# Original code by Michael Christoffersen, Brandon Tober
# Last edited 2019-11-22 by Stefano Nerozzi (stefano.nerozzi@utexas.edu)

import numpy as np
import sys
import math
import warnings

# Parameters
marsradius = 3396  # Mars radius in km
nd_value = 9999.999 # no data value
csize_in = 5 # size of input grid cells in degrees
csize_out = 2 # size of output grid cells in degrees
avg_dist = 2.0  # search radius [sigma]
halfd = 265 # 50% signal radius [km]
infile = 'cl_5x5.csv'
outfile = 'cl_5x5_gauss_2s.csv'
in_dir = '/mnt/d/Dropbox/MARS/GRS/elements-v1/unsmoothed/5x5/'
out_dir = '/mnt/d/Dropbox/MARS/GRS/elements-v1/processed2/5x5/'
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
  data = np.genfromtxt(infile, delimiter=',')

  # replace no data value with NaN
  data[data==nd_value] = np.NaN

  # Put input data into a grid
  latcells_in = int(180/csize_in) # number of input latitude cells
  loncells_in = int(360/csize_in) # number of input longitude cells
  datagrid = np.zeros((latcells_in, loncells_in, 3))
  datagrid[:] = np.NaN
  for k in range(len(data)):
    i = int( (90 - data[k,0] - csize_in/2) / csize_in )
    j = int( (data[k,1] - csize_in/2 ) / csize_in )
    datagrid[i,j,0:] = data[k,2:]

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
  lonl = np.matrix(np.linspace(csize_in/2, 360-csize_in/2, loncells_in))
  latl = np.matrix(np.linspace(90-csize_in/2, -90+csize_in/2, latcells_in))
  lon = np.dot(np.ones((datagrid.shape[0],1)), lonl)
  lat = np.dot(np.transpose(latl),np.ones((1,datagrid.shape[1])))

  # Make n*n degree output grid
  latcells_out = int(180/csize_out) # number of output latitude cells
  loncells_out = int(360/csize_out) # number of output longitude cells
  outgrid = np.zeros((latcells_out, loncells_out, 3))
  olon = np.linspace(csize_out/2, 360-csize_out/2, loncells_out)
  olat = np.linspace(90-csize_out/2, -90+csize_out/2, latcells_out)

  # Do averaging, iterating over output grid
  for i in range(latcells_out):
    for j in range(loncells_out):
      # Haversine distance from each point to rest of grid
      dist = havs(olat[i], olon[j], lat, lon, marsradius)
      # Calculate weights for average based on gaussian function
      weights = math.e ** ((-dist**2) / (2*sigma**2))
      # Mask out points and weights outside search radius
      subpts = datagrid[dist < avg_dist*sigma]
      weights = np.transpose([weights[dist < avg_dist*sigma],]*3)
      # Create masked arrays to avoid NaNs later
      subpts = np.ma.masked_array(subpts, np.isnan(subpts))
      weights = np.ma.masked_array(weights, np.isnan(subpts))
      # Calculate average ignoring NaNs
      outgrid[i,j] = np.ma.filled(np.ma.average(subpts, axis=0, weights=weights), np.nan)
      # Propagate uncertainty (averaging improves uncertainty by sqrt(number of averaged samples))
      outgrid[i,j,1:] = np.divide( outgrid[i,j,1:], np.sqrt(np.sum(weights[:,1:], axis=0)))
      # Print progress
      progress = int(i/latcells_out*100)
      sys.stdout.write('\r')
      sys.stdout.write("[%-50s] %d%%" % ('='*(progress//2), progress))
      sys.stdout.flush()
  
  # Go to new line after printing progress bar
  print("\r\n")

  # Write out text file in the same format as the original, without NaN/no data values
  with open(outfile, "w") as f:
    for i in range(latcells_out):
      for j in range(loncells_out):
        if ~np.isnan(outgrid[i,j,0]):
          f.write("{},{},{:.8f},{:.8f},{:.8f}\n".format(olat[i], olon[j], outgrid[i,j,0], outgrid[i,j,1], outgrid[i,j,2], ))
  print("Done! Output is in: " + out_dir)
   
main()