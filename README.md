# Nanopore-TnSeq-Pipeline
This repository contains a custom pipeline developed for the analysis of transposon sequencing (Tn-seq) data generated using Oxford Nanopore Technologies. 

## Pipeline Overview

1. **Read filtering**: Retain reads 150–180 bp in length.
2. **Transposon screening**: Select reads containing a forward or reverse transposon adapter sequence.
3. **Trimming**: Extract the 14 bp region immediately downstream of the adapter.
4. **Alignment**: Map trimmed reads to a reference genome using Bowtie2. Only unique, perfect matches are retained.
5. **TA-site quantification**: Count insertions at each TA site and export to `.wig` format.
6. **WIG harmonization**: Standardize `.wig` files by adding missing positions (with 0 counts) for consistent comparison.

## Scripts

This pipeline includes three core scripts:

- **`tnseq_prepmap.py`**:  
  Filters reads by length and adapter presence, trims them, and aligns insert sequences to a reference genome using Bowtie2. Outputs `.wig` files and sample-specific statistics.

- **`tnseq_combine_stats.py`**:  
  Combines multiple stats files (produced by `tnseq_prepmap.py`) into a single summary table for easy comparison across samples.

- **`tnseq_harmonize_wigs.py`**:  
  Ensures all `.wig` files share the same genomic positions by adding missing sites with zero counts—essential for downstream comparisons and visualization.


## Example Usage

```bash
python tnseq_prepmap.py -i reads.fastq -o sample1 -r reference.fasta
python tnseq_combine_stats.py *_stats.txt
python tnseq_harmonize_wigs.py *.wig
```

## Citation

If you use this pipeline in your research, please cite the following link:
https://github.com/MariaJuanpereBorras/Nanopore-TnSeq-Pipeline

## License

This project is licensed under the [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) license.

