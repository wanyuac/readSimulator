#!/usr/bin/env python

"""
This script simulates paired-end short reads from reference sequences. It is inspired by Ryan's script
make_synthetic_reads.py. At present, it supports simulators wgsim and ART. A user can simulate reads based
on circular topology (via setting --iterations >1) of a completed bacterial genome or on linear topology
(setting --iterations 1) of contigs. I think it is safer to run this script than using wgsim directly for a
multi-FASTA file because I am not sure whether wgsim will treat sequences as independent or concatenate them.

Arguments:
The input reference FASTA file may contain multiple contigs (that is, a multi-FASTA file). Reads simulated
from contigs of the same FASTA file will be pooled into the same paired-end read set.
'opts' is a string that will be passed to the simulator as options.
'iterations'
    By default, it equals one, which indicates the script to treat every reference sequence as linear.
    Reference sequences are treated as circular when this argument is greater than one (say, n). Herein,
    the script rotates every reference sequence conceptually through changing its start position and
    concatenating sequences beside this arteficial start point. The simulator generates reads from this
    recombined sequence. This process is repeated for n times so as to get an even coverage across the
    whole sequence.
Simulator must be either 'wgsim' or 'art'.
check_cmd: configure this option to print command lines without actual execution.

Example commands:
    python readSimulator.py --input ../references/*.fasta --simulator wgsim --simulator_path ../apps/wgsim \
    --outdir shredded_reads --iterations 10 --readlen 76 --opts '-e 0 -r 0 -R 0 -X 0 -h -S 5' > readsim.log
    
    python readSimulator.py --input ../references/*.fasta --simulator art \
    --simulator_path ../apps/art/art_illumina \
    --outdir shredded_reads --iterations 10 --readlen 100 \
    -opts '-m 300 -s 25 -qU 35 -nf 0 -rs 10 -ss HS20 -qs -1.5 -qs2 -2 -na' --check_cmd > readsim.log

For some cluster users, wgsim may be compiled as a component of SamTools. In this case, one can simply call
wgsim after loading the environment module SAMtools.

Limitations:
1. This script does not support variable relative depths for each contig because the number of contigs varies
across input files.
2. It does not support simulation of PacBio reads yet.

Python version: 3.5.2+
Author: Yu Wan (wanyuac@126.com)
Licence: GNU GPL 2.1
Publication: 23-24/9/2016
Last modification: 14/5/2021
"""

from argparse import ArgumentParser
import os, sys, random, subprocess

SUPPORTED_SIMULATORS = ['wgsim', 'art']

def parse_arguments():
    
    parser = ArgumentParser(description = "Simulate reads based on reference sequences")
    parser.add_argument('--input', nargs = '+', type = str, required = True,\
                        help = "FASTA files containing reference sequences for shredding")
    parser.add_argument('--simulator', type = str, default = 'wgsim', required = False,\
                        help = "Name for the simulator. Supported: %s." % ', '.join(SUPPORTED_SIMULATORS))
    parser.add_argument('--simulator_path', type = str, default = 'wgsim', required = False,\
                        help = "The path and file name of the shredder")
    parser.add_argument('--outdir', type = str, default = '.', required = False,\
                        help = "Directory for output files")
    parser.add_argument('--iterations', type = int, default = 1, required = False,\
                        help = "Number of arbitrary sequence breakpoints for read simulation. (default: 1, treating the template genome as a linear sequence)")
    parser.add_argument('--readlen', type = int, default = '100', required = False,\
                        help = "Read length")
    parser.add_argument('--depth', type = float, default = '70', required = False,\
                        help = "Required approximate fold coverage of every sequence. Default: 70 folds")
    parser.add_argument('--opts', type = str, required = True, \
                        help = "A string of other options and arguments to be passed to the read simulator.")
    parser.add_argument('--check_cmd', action = 'store_true', required = False,\
                        help = "Configure it to print command lines only.")
    
    return parser.parse_args()

def main():
    
    args = parse_arguments()
    run = not args.check_cmd
    
    # check if the simulator is set correctly
    tool = args.simulator.lower()
    if tool in SUPPORTED_SIMULATORS:
        simulator = {tool : args.simulator_path}  # make a dictionary for the convenience to call the simulator
        print("Simulator: " + tool + " at " + args.simulator_path)
    else:
        sys.exit("Error: the simulator must be selected from " + "/".join(SUPPORTED_SIMULATORS) + "!")
    
    # make an output folder
    if run:
        if not os.path.exists(args.outdir):
            os.system('mkdir ' + args.outdir)
        
    # loop through every sample
    for fasta in args.input:
        sample = os.path.splitext(os.path.basename(fasta))[0]  # extract sample name from the file name
        seq_names, ref_seqs = load_fasta(fasta)  # return a list of sequence names and a list of sequences in the current file
        
        # initialise output files
        readsets = [os.path.join(args.outdir, sample + '_' + '1.fastq'),\
                    os.path.join(args.outdir, sample + '_' + '2.fastq')]
        if run:
            for _, f in enumerate(readsets):
                open(f, 'w').close()  # create new files or erase existing content of output files
        
        # simulate reads for every contig in the current file
        for j, seq in enumerate(ref_seqs):
            if args.iterations > 1:  # circular mode
                for i in list(range(0, args.iterations)):  # change arteficial breakpoints on every circular molecule
                    simulate_reads(sample = sample, seq_name = seq_names[j], seq = seq,\
                                   readlen = args.readlen, depth = args.depth / args.iterations,\
                                   filenames = readsets, simulator = simulator, opts = args.opts,\
                                   random_start = True, run = run)
            else:  # linear mode
                simulate_reads(sample = sample, seq_name = seq_names[j], seq = seq,\
                               readlen = args.readlen, depth = args.depth,\
                               filenames = readsets, simulator = simulator, opts = args.opts,\
                               random_start = False, run = run)
        
        # finally, compress read sets of the current sample
        if run:
            for _, f in enumerate(readsets):
                zip_file(f)
            print("Simulation for sample %s is finished successfully. Simulated reads are stored in %s and %s" %\
                  (sample, readsets[0] + ".gz", readsets[1] + ".gz"))

def simulate_reads(sample, seq_name, seq, readlen, depth, filenames, simulator, opts, random_start, run):
    
    # make a temporal file as the reference for read simulation
    if random_start:  # Circular mode: randomly rotate the sequence and save it to a temporary file.
        random_start = random.randint(0, len(seq) - 1)
        refseq = seq[random_start : ] + seq[ : random_start]  # rotate the sequence by breaking the circle at random_start so that a new sequence starts here
    else:
        refseq = seq  # linear mode
    refseq_tmp = sample + '_' + 'refseq.fasta'  # add a prefix to the common temporary file names so that multiple jobs can be executed simultaneously
    if run:
        refseq_file = open(refseq_tmp, 'w')  # create a new temporary FASTA file
        refseq_file.write('>' + seq_name + '\n')
        refseq_file.write(refseq + '\n')
        refseq_file.close()  # Nothing gets writen into this file unless it is closed. Otherwise, the reference sequence for wgsim is null.

    # simulate short reads with wgsim
    read_pair_num = int(round(depth * len(seq) / (2 * readlen), 0))  # for paired-end reads
    tool = list(simulator.keys())[0]
    read1_tmp = sample + '_tmp_1.fq'  # ART only produces files with a suffix of "fq". Wgsim accepts any suffices.
    read2_tmp = sample + '_tmp_2.fq'
    
    if tool == 'wgsim':
        cmd = [simulator[tool], '-N', str(read_pair_num), '-1', str(readlen), '-2', str(readlen)] + \
        opts.split(" ") + [refseq_tmp, read1_tmp, read2_tmp]
    elif tool == 'art':
        cmd = [simulator[tool], '--in', refseq_tmp, '--rcount', str(read_pair_num), \
               '--len', str(readlen), '--paired', '--out',  sample + '_tmp_'] + opts.split()
    else:  # support to other tools
        pass
    
    print(" ".join(cmd))  # print the command line to stdout
    if run:
        process = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        out, err = process.communicate()
    
        # append new synthetic reads to output files before compression
        readset_1 = open(filenames[0], "a")
        readset_2 = open(filenames[1], "a")
    
        with open(read1_tmp, "r") as f1:
            for line in f1:
                readset_1.write(line) 
        with open(read2_tmp, "r") as f2:
            for line in f2:
                readset_2.write(line)
        
        # delete temporal files
        os.remove(read1_tmp)
        os.remove(read2_tmp)
        os.remove(refseq_tmp)

def load_fasta(filename): # type: (str) -> list[tuple[str, str]]
    # Returns paired sequence names and sequences for each fasta file.
    
    seq_names = []  # a list of sequence IDs
    seqs = []
    fasta_file = open(filename, "rU")
    name = ""
    sequence = ""
    
    with open(filename, "rU") as f:
        content = f.read().splitlines()  # Empty lines ('\n') will not go into the list
    
    for line in content:
        if line.startswith('>'): # Each header line leads information of a new circular replicon.
            if name:
                # store previous sequence information when the pointer encounters a new sequence
                seq_names.append(name.split(" ")[0])  # parse last FASTA definition line and only keep the sequence ID
                seqs.append(sequence)  # add last sequence into the list
                sequence = ""
            name = line[1:]  # store a new FASTA definition line without the first character ('>')
        else:
            sequence += line  # concatenate a new line with previous sequences
    if name:  # store the last name and sequence
        seq_names.append(name.split(" ")[0])
        seqs.append(sequence)
        
    return seq_names, seqs  # two parallel lists

def zip_file(filename):
    command = ["gzip", filename]
    process = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    _, _ = process.communicate()

if __name__ == '__main__':
    main()
