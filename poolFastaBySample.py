"""
This script pools FASTA files of the same sample together to make a large multi-FASTA file.
It can be used for preparing input files of readSimulator_slurm.py (multi-sample) and readSimulator.py (single sample).
It expects input FASTA files to follow the filename convension used in my download_NCBI_records.py (https://github.com/wanyuac/BINF_toolkit):
field1(delimiter)field2.extension. For instance, strain1__LFJS01000001.fasta is a valid filename, which is comprised of a strain name
and an NCBI accession.

An input file must consist of two columns (tab-delimited) without column names. The first column may contain sample/strain/isolate
names and the second column is comprised of identifiers (comma-delimited) for the field 2.
For instance, a row in an input file reads:
    NJST258_1\tCP006923,CP006924,CP006927,CP006928,CP006925,CP006926
Accordingly, NJST258_1__CP006923.fasta ... NJST258_1__CP006927.fasta will be pool together into NJST258_1.fasta.

Arguments
--input: a text file following the aforementioned structure.
--filename_delim: the delimiter for two fields in input file names (default: "__").
--col_delim: the delimiter separating two columns of the input file (default: "\t").
--col2_delim: the delimiter for values in the second column of the input file (default: ",").
--indir: the directory where input FASTA files are stored.
--input_ext: filename extension of input FASTA files (default: fasta)
--outdir: the output directory
--output_ext: filename extension of output FASTA files (default: fasta)

Example:
python poolFastaBySample.py --input strain_accessions.txt --filename_delim "__" --col_delim "\t" --col2_delim "," \
indir genbank/fasta --input_ext fasta --outdir genbank/fasta_pooled --output_ext fna

Python version: 3.5.2
Author: Yu Wan (wanyuac@gmail.com)
Licence: GNU GPL 2.1
Development history: 14 April 2017
"""

from argparse import ArgumentParser
import os, sys

def parse_arguments():
    parser = ArgumentParser(description = "Pool FASTA files together by each sample.")
    parser.add_argument("--input", type = str, required = True, help = "A two-column text file for filenames")
    parser.add_argument("--filename_delim", type = str, default = "__", required = False, help = "Delimiter for two fields in input file names")
    parser.add_argument("--col_delim", type = str, default = "\t", required = False, help = "Delimiter separating two columns of the input file")
    parser.add_argument("--col2_delim", type = str, default = ",", required = False, help = "Delimiter for values in the second column of the input file")
    parser.add_argument("--indir", type = str, default = ".", required = False, help = "Directory of input FASTA files")
    parser.add_argument("--input_ext", type = str, default = "fasta", required = False, help = "Filename extension of input FASTA files")
    parser.add_argument("--outdir", type = str, default = "fasta_pooled", required = False, help = "Output directory")
    parser.add_argument("--output_ext", type = str, default = "fasta", required = False, help = "Filename extension of output FASTA files")
    
    return parser.parse_args()

def parse_names(f, col_delim, col2_delim):
    # parse the two-column input file
    with open(f, "rU") as name_table:
        lines = name_table.read().splitlines()
    names = {}
    for line in lines:
        key, vals = line.split(col_delim)  # separate column 1 and column 2
        names[key] = vals.split(col2_delim)  # {key : [id 1, id 2, ..., id n]}
    
    return(names)

def main():
    args = parse_arguments()
    
    if not os.path.exists(args.outdir):  # check and set up the output directory
        os.system("mkdir " + args.outdir)
    
    fastas = parse_names(f = args.input, col_delim = args.col_delim, col2_delim = args.col2_delim)
    for sample, accessions in fastas.items():
        target_files = [args.filename_delim.join([sample, acc]) for acc in accessions]  # sample__accession1, ..., sample_accessionN
        target_files = [".".join([name, args.input_ext]) for name in target_files]
        target_files = [os.path.join(args.indir, name) for name in target_files]
        os.system(" ".join(["cat", " ".join(target_files), ">", os.path.join(args.outdir, sample + "." + args.output_ext)]))
    
if __name__ == '__main__':
    main()
    