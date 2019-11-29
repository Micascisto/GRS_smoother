# GRS Smoother

## Description

Simple set of python scripts to smooth Mars PDS Gamma Ray Spectrometer data
Original code by Brandon Tober and Michael Christoffersen

## Some useful info

- *readDat.py* parses .TAB files from PDS and creates much friendlier csv files.

- *withinRad.py* is ~10x slower than *rollmean.py*, and produces a slightly shifted output

- Structure of latitude and longitude arrays used by *roolmean.py*:
```
    lat 5x5 deg array looks like
    [ 87.5, 87.5, 87.5, ... ]
    [ 82.5, 82.5, 82.5, ... ]
    [ 77.5, 77.5, 77.5, ... ]
    [ ...., ...., ...., ... ]

    lon 5x5 deg array looks like
    [ 2.5, 7.5, 12.5, ...]
    [ 2.5, 7.5, 12.5, ...]
    [ 2.5, 7.5, 12.5, ...]
    [ ..., ..., ...., ...]
```
## Notes, future improvements

- *rollmean.py* is pretty fast, but would probably be faster if parallelized
- Correct Mars radius for polar flattening when calculating distances
