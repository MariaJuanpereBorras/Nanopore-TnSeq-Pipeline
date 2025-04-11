#!/usr/bin/env python3
# todoc

"""
tnseq_harmonize_wigs
---------------
.. module:: tnseq_harmonize_wigs.py
  :synopsis: add missing genomics locations to wig files
.. moduleauthor:: Maria Juanpere Borras and Jos Boekhorst

Script to add missing genomics locations to wig files, ensuring they have exactly the same sites (as required by Transit track viewer)

Typical run::

    tnseq_harmonize_wigs.py *.wig

Run the script with '-h' for a list of options.

"""
import glob
import argparse
import sys


def read_wig(filename):
    data = {}
    with open(filename, 'r') as f:
        lines = f.read().strip()
    lines = lines.split('\n')
    data['header'] = lines[0]
    data['counts'] = {}
    locations = set([])
    for line in lines[1:]:
        lineg = line.split()
        location = int(lineg[0])
        locations.add(location)
        count = int(lineg[1])
        data['counts'][location] = count
    return data, locations


def main(options):
    infiles = options['other_arguments']
    # guess contig names from filenames; assumes no "_" in contig names, which won't work for many NCBI genomes!
    contigs = set([])
    for name in infiles:
        contig = name.split('/')[-1].split('_')[-1].replace('.wig', '')
        contigs.add(contig)
    contigs = list(contigs)
    print(f"Found contigs: {'; '.join(contigs)}")

    for contig in contigs:
        wigs = {}
        infiles = glob.glob(f"*{contig}.wig")
        all_locations = set([])
        for name in infiles:
            wigs[name], locations = read_wig(name)
            all_locations = all_locations | locations
        all_locations = list(all_locations)
        all_locations.sort()
        for name in infiles:
            with open(name.replace('.wig', '_harmonized.wig'), 'w') as w:
                w.write(wigs[name]['header'] + '\n')
                for location in all_locations:
                    try:
                        count = wigs[name]['counts'][location]
                    except KeyError:
                        count = 0
                    w.write(f"{location}  {count}\n")
    print("Done")



############
# SETTINGS #
############
# note: default parameters are set in argparse object (top of __main__)
description = "Script to add missing genomics locations to wig files, ensuring they have exactly the same sites (as required by Transit track viewer)"

# main program
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description, add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('other_arguments', help="individual unnamed parameters", nargs='+')
    options = vars(parser.parse_args())
    main(options)




    
