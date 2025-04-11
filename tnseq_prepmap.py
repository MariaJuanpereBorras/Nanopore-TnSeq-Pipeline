#!/usr/bin/env python3
# todoc

"""
tnseq_prepmap
---------------
.. module:: tnseq_prepmap
  :synopsis: Filter on length, trim adapter and map to reference for tnseq analysis
.. moduleauthor:: Maria Juanpere Borras and Jos Boekhorst

Filter and trim fastq file upto as specific sequence. Nothing fancy, does NOT allow any mismatches etc.

Typical run::

    tnseq_prepmap.py -i my_reads.fastq -o my_reads -a ACTTATCATCCAACCTGTTA

Run the script with '-h' for a list of options.

"""

import argparse
import sys
import subprocess
from Bio import SeqIO


def write_stats(options, stats):
    with open(options['prefix'] + '_stats.txt', 'w') as f:
        f.write(f"command: {' '.join(sys.argv)}\n")
        for element in stats:
            f.write(f"{element}: {stats[element]}\n")


def do_command(command):
    sys.stderr.write(f"Executing: {command}\n")
    returned_value = subprocess.call(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    if returned_value != 0:
        sys.stderr.write("Command returned non-zero exit status, aborting\n")
        sys.exit(1)


def find_adapter(sequence, adapter, stats, mismatches=0):

    if mismatches == 0:
        adapter_start = sequence.find(adapter)
        if adapter_start !=-1:
            stats['with_perfect_adapter'] += 1
        return adapter_start, stats

    if mismatches == 1:
        adapter_start = sequence.find(adapter)
        if adapter_start != -1: # if already a no-mismatch hit, return that one
            stats['with_perfect_adapter'] += 1
            return adapter_start, stats
        nts = "AGCT"
        all_variants = set([])
        for position in range(len(adapter)):
            for nt in nts:
                variant = adapter[:position] + nt + adapter[position+1:]
                all_variants.add(variant)
        all_variants = list(all_variants)
        all_variants.sort()
        for variant in all_variants:
            pos = sequence.find(variant)
            if pos != -1:
                stats['with_mismatched_adapter'] += 1
                return pos, stats
        return -1, stats

    sys.stderr.write("Only 0 or 1 mismatches currently supported\n")
    sys.exit(1)


def main(options):

    stats = {'total_reads': 0,
             'correct_length': 0,
             'with_perfect_adapter': 0,
             'with_mismatched_adapter': 0,
             'total_mapped': 0,
             'perfect_matches': 0,
             'unique_perfect_matches': 0}

    sys.stderr.write(f"Reading {options['fastq']}\n")
    fastq = SeqIO.parse(options['fastq'], format="fastq")

    # lenght & adapter fitlering
    sys.stderr.write("Filtering and adapter trimming\n")
    insert_outname = options['prefix'] + '_inserts.fastq'
    with open(insert_outname, 'w') as f:
        for element in fastq:
            stats['total_reads'] += 1
            if int(options['minlen']) <= len(element.seq) <= int(options['maxlen']):
                stats['correct_length'] += 1
                adapter_start, stats = find_adapter(element.seq, options['adapter'], stats, mismatches=options['mismatches'])
                if adapter_start == -1:  # not found; either missing or reverse complement
                    tmp_id = element.id + "_RC"  # ID is lost when reverse complementing
                    tmp_description = element.description
                    element = element.reverse_complement()
                    element.id = tmp_id
                    element.description = tmp_description
                    adapter_start, stats = find_adapter(element.seq, options['adapter'], stats, mismatches=options['mismatches'])
                if adapter_start != -1:  # found
                    insert_start = adapter_start + len(options['adapter'])
                    insert_stop = insert_start + int(options['insertsize'])
                    element = element[insert_start:insert_stop]
                    SeqIO.write(element, f, "fastq")

    # bowtie2 stuff
    bowtie_outname = options['prefix'] + "_mapping.sam"
    do_command(f"bowtie2 -x {options['fna']} -U {insert_outname} -S {bowtie_outname}")

    # parsing bowtie2 output (=SAM with bowtie2 tags)
    coverage = {}
    with open(bowtie_outname, 'r') as f:
        for line in f:
            if line[0] != "@":  # skip the comment lines
                lineg = line.rstrip('\n').split('\t')
                # query = lineg[0]
                flag = lineg[1]
                if flag != "4":  # 4: no reported alginments
                    stats['total_mapped'] += 1
                    target = lineg[2]
                    position = int(lineg[3])
                    # quality = int(lineg[4])
                    tags = lineg[11:]
                    if target not in coverage:
                        coverage[target] = {}

                    if 'XM:i:0' in tags:  # no mismatches
                        stats['perfect_matches'] += 1
                        if not ('AS:i:0' in tags and 'XS:i:0' in tags):  # no other perfect hit
                            stats['unique_perfect_matches'] += 1
                            if position not in coverage[target]:
                                coverage[target][position] = 0
                            coverage[target][position] += 1

    #  write to wig files
    for target in coverage:
        with open(f"{options['prefix']}_{target}.wig", 'w') as f:
            positions = list(coverage[target].keys())
            positions.sort()
            f.write(f"variableStep chrom={target}\n")
            for position in positions:
                f.write(f"{position}  {coverage[target][position]}\n")

    # conver gff to prot_table
    # do_command(f"transit convert gff_to_prot_table {options['gff']} {options['prefix']}.prot_table")

    stats['with_adapter'] = stats['with_perfect_adapter'] + stats['with_mismatched_adapter']
    write_stats(options, stats)
    sys.exit()


############
# SETTINGS #
############
# note: default parameters are set in argparse object (top of __main__)
description = "Filter on length, trim adapter and map to reference for tnseq analysis"

# main program
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=description, add_help=True,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', dest='fastq', help='input fastq file', required=True)
    parser.add_argument('-o', dest='prefix', help='prefix for output file', required=True)
    parser.add_argument('-r', dest='fna', help='reference genome (fasta) file', required=True)
    parser.add_argument('-a', dest='adapter', help='trim upto this sequence', default="ACTTATCATCCAACCTGTTA")
    parser.add_argument('-m', dest='minlen', help='minimum sequence length', default=150)
    parser.add_argument('-n', dest='maxlen', help='maximum sequence length', default=180)
    parser.add_argument('-s', dest='insertsize', help='insert size', default=14)
    parser.add_argument('-x', dest='mismatches', help='adapter mismatches allowed', default=0, type=int)
    options = vars(parser.parse_args())
    main(options)
