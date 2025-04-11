#!/usr/bin/env python3
# todoc

"""
tnseq_combine_stats
---------------
.. module:: tnseq_combine_stats
  :synopsis: Combine tnseq stats file
.. moduleauthor:: Maria Juanpere Borras and Jos Boekhorst

ombine tnseq stats file

Typical run::

    tnseq_combine_stats.py *_stats.txt

Run the script with '-h' for a list of options.

"""

import argparse
import sys

############
# SETTINGS #
############
# note: default parameters are set in argparse object (top of __main__)
description = "ombine tnseq stats file"


# main program
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description, add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('other_arguments', help="individual unnamed parameters", nargs='+')
    options = vars(parser.parse_args())

    data = {}

    for name in options['other_arguments']:
        data[name] = {}
        with open(name) as f:
            labels = []
            for line in f:
                line = line.strip()
                label, value = line.split(': ')
                data[name][label] = value
                labels.append(label)

    print('\t'.join(['sample'] + labels))
    for name in options['other_arguments']:
        line = [name]
        for label in labels:
            line.append(data[name][label])
        print('\t'.join(line))

    sys.exit()
