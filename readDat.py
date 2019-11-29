# Read .tab data file and remove white space to export as csv

import re, os, argparse

def get_args():
    """Function to parse arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", type=str, required=True, help="Input directory with absolute path")
    args = parser.parse_args()
    directory = args.input

    return directory

# Runs through files under directory
directory = get_args()
for infile in os.listdir(directory):
    if infile.endswith('tab'):
        infile = os.path.join(directory, infile)

        # Open file and read all lines
        with open(infile) as f:
            dat = f.readlines()

        # Empty list to hold modified data
        dat_csv = []

        # For each line, remove whitespace and replace with commas for csv format - use re (regular expression)
        for line in dat:
            lines = re.sub("\s+", ",", line.strip())
            # Append modified line to new list
            dat_csv.append(lines)

        # Output data to csv with newline character following each line
        outfile = infile.rstrip('.tab') + '.csv'
        with open(outfile,'w') as f:
            for line in dat_csv:
                f.write(line)
                f.write('\n')