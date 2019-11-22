import re, os
"""
read .tab data file and remove white space to export as csv
"""
directory = '/mnt/d/Dropbox/MARS/GRS/elements-v1/unsmoothed/2x2/'

# runs through files under directory
for infile in os.listdir(directory):
    if infile.endswith('tab'):
        infile = directory+infile

        # open file and read all lines
        with open(infile) as f:
            dat = f.readlines()

        # empty list to hold modified data
        dat_csv = []

        # for each line, remove whitespace and replace with commas for csv format - use re (regular expression)
        for line in dat:
            lines = re.sub("\s+", ",", line.strip())
            # append modified line to new list
            dat_csv.append(lines)

        # output data to csv with newline character following each line
        outfile = infile.rstrip('.tab') + '.csv'
        with open(outfile,'w') as file:
            for line in dat_csv:
                file.write(line)
                file.write('\n')