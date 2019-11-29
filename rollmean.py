# Pyhton script to smooth PDS Gamma Ray Spectrometer grids
# Original code by Michael Christoffersen, Brandon Tober
# Last edited 2019-11-22 by Stefano Nerozzi (stefano.nerozzi@utexas.edu)

# Usage: rollmean.py -i <input file> -o <output file> -a <averaging radius> -c <output cell size>

import sys, os.path, math, argparse
import numpy as np

# Constants
MARS_RADIUS = 3389.5  # Mars radius in km
ND_VALUE = 9999.999 # "No data" value
ELEMENTS = {
    "si": 219,
    "fe": 267,
    "h": 222,
    "cl": 270,
    "k": 215,
    "th": 240,
}

def get_args():
    """Function to parse arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="Input file with absolute path")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output file with absolute path")
    parser.add_argument("-a", "--average-radius", type=float, required=True, help="Averaging radius (number of standard deviations)")
    parser.add_argument("-c", "--output-cell-size", type=float, required=True, help="Output cell size")
    args = parser.parse_args()
    infile = args.input
    outfile = args.output
    avg_dist = args.average_radius  # search radius [sigma]
    csize_out = args.output_cell_size # size of output grid cells in degrees

    return infile, outfile, avg_dist, csize_out

# Haversine formula for distance
def havs(latp, lonp, lats, lons):
    # latp, lonp are scalar
    # lats, lons are 2D arrays
    # Output is array with distance from (latp, lonp) to every point in (lats, lons)
    lats = np.array(lats)
    lons = np.array(lons)
    dLat = np.radians(np.subtract(lats, latp))
    dLon = np.radians(np.subtract(lons, lonp))
    latp = np.radians(latp)
    lats = np.radians(lats)
 
    a = np.power(np.sin(dLat/2),2) + np.cos(latp) * np.cos(lats) * np.power(np.sin(dLon/2),2)
    c = 2 * np.arcsin(np.sqrt(a))
 
    return np.multiply(MARS_RADIUS, c)

def main():
    # Get arguments
    infile, outfile, avg_dist, csize_out = get_args()

    # Process some stuff
    filename = os.path.basename(infile)
    element = filename.split("_")[0] # read element symbol
    csize_in = int(filename.split("_")[1].split("x")[0]) # size of input grid cells in degrees
    hwhm = ELEMENTS[element] # find 50% signal radius [km] (half width half maximum)
    sigma = hwhm / np.sqrt(2*np.log(2)) # Calculate 1 sigma radius; HWHM = sigma * sqrt(2*ln(2))

    # Print some useful info
    print("Input file: {}".format(infile))
    print("50% signal radius: {} km; 1-sigma radius: {:.0f} km".format(hwhm, sigma))
    print("Averaging radius: {:.0f} km or {} sigma".format(avg_dist*sigma, avg_dist))

    # Load input data file, if it exists
    if not os.path.isfile(infile):
        print("Fucking bananas!!!")
        sys.exit()
    data = np.genfromtxt(infile, delimiter=",")

    # replace no data value with NaN
    data[data==ND_VALUE] = np.NaN

    # Put input data into a grid
    latcells_in = int(180/csize_in) # number of input latitude cells
    loncells_in = int(360/csize_in) # number of input longitude cells
    datagrid = np.zeros((latcells_in, loncells_in, 3))
    datagrid[:] = np.NaN
    for row in data:
        i = int( (90 - row[0] - csize_in/2) / csize_in )
        j = int( (row[1] - csize_in/2 ) / csize_in )
        datagrid[i,j,0:] = row[2:]

    # Generate lat and lon arrays, used for haversine distance calculation
    lonl = np.matrix(np.linspace(csize_in/2, 360-csize_in/2, loncells_in))
    latl = np.matrix(np.linspace(90-csize_in/2, -90+csize_in/2, latcells_in))
    lon = np.dot(np.ones((datagrid.shape[0],1)), lonl)
    lat = np.dot(np.transpose(latl),np.ones((1,datagrid.shape[1])))

    # Make output grid
    latcells_out = int(180/csize_out) # number of output latitude cells
    loncells_out = int(360/csize_out) # number of output longitude cells
    outgrid = np.zeros((latcells_out, loncells_out, 3))
    olon = np.linspace(csize_out/2, 360-csize_out/2, loncells_out)
    olat = np.linspace(90-csize_out/2, -90+csize_out/2, latcells_out)

    # Calculate weighted average, iterating over output grid
    for i in range(latcells_out):
        for j in range(loncells_out):
            # Haversine distance from each point to rest of grid
            dist = havs(olat[i], olon[j], lat, lon)
            # Calculate weights for average based on gaussian function
            weights = np.power(math.e, -0.5*np.power(dist/sigma, 2))
            # Mask out points and weights outside search radius
            distmask = dist < avg_dist*sigma
            subpts = datagrid[distmask]
            weights = np.transpose([weights[distmask],]*3)
            # Create masked arrays to avoid NaNs later
            nanmask = np.isnan(subpts)
            subpts = np.ma.masked_array(subpts, nanmask)
            weights = np.ma.masked_array(weights, nanmask)
            # Calculate average ignoring NaNs
            outgrid[i,j] = np.ma.filled(np.ma.average(subpts, axis=0, weights=weights), np.nan)
            # Propagate uncertainty (averaging improves uncertainty by sqrt(number of averaged samples))
            outgrid[i,j,1:] = np.divide(outgrid[i,j,1:], np.sqrt(np.sum(weights[:,1:], axis=0)))
            # Print progress
            progress = int(i/latcells_out*100+1)
            sys.stdout.write("\r")
            sys.stdout.write("[%-50s] %d%%" % ("="*(progress//2), progress))
            sys.stdout.flush()
    
    # Go to new line after printing progress bar
    print("\r")

    # Write out text file in the same format as the original, without NaN/no data values
    with open(outfile, "w") as f:
        for i in range(latcells_out):
            for j in range(loncells_out):
                if ~np.isnan(outgrid[i,j,0]):
                    f.write("{},{},{:.8f},{:.8f},{:.8f}\n".format(olat[i], olon[j], outgrid[i,j,0], outgrid[i,j,1], outgrid[i,j,2], ))
    print("Done! Output is in: " + outfile)

main()