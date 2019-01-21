"""
This script submits SLURM jobs of readSimulator.py for multiple reference sequences so as to save time.
At present, it only supports the VLSCI cluster Helix, but it is easy to adapt this script for other SLURM-based systems.

Arguments:
    --input: reference sequences in the FASTA format for simulation
    --opts: a string of options and arguments that will be pasted directly into the command line of readSimulator.py.
    --path: the directory where the readSimulator.py locates. By default, it is the current working directory.
    --mem: memoery (Mb) assigned to each job. Default: 2048 Mb.
    --wait: seconds before submitting the next job. Increasing the interval to avoid flooding the queuing system with
    too many jobs, which may result in a loss of outputs.
    --walltime: the maximal run time of each job. Default: 0-0:30:0 (30 min).

Example:
    python readSimulator_slurm.py --input *.fasta --path codes --mem 4096 --walltime 600 \
    --opts "--simulator wgsim --simulator_path wgsim --outdir shredded_reads --iterations 10 \
    --readlen 100 --depth 70 --opts '-e 0 -r 0 -R 0 -X 0 -S 5 -h'"

Python version: 3.5.2
Author: Yu Wan (wanyuac@gmail.com)
Licence: GNU GPL 2.1
First edition: 24 Sep 2016; the latest edition: 20 May 2018
"""

from argparse import ArgumentParser
import os
import sys
import time

# Cluster specific arguments. Change them to adapt to another SLURM system. ###############
# Helix
PARTITION = 'sysgen,main'
MODULE_PYTHON = 'Python/3.5.2-vlsci_intel-2015.08.25'
MODULE_SAMTOOLS = 'SAMtools/1.3.1-vlsci_intel-2015.08.25-HTSlib-1.3.1'

# Functions ###############

def parse_arguments():
    parser = ArgumentParser(description = "Simulate reads based on reference sequences")
    parser.add_argument('--input', nargs = '+', type = str, required = True,\
                        help = "FASTA files containing reference sequences for shredding")
    parser.add_argument('--path', type = str, default = '.', required = False,\
                        help = "The directory containing the readSimulator.py script")
    parser.add_argument('--mem', type = str, default = '2048', required = False,\
                        help = "Memoery (Mb) per job. Default: 2048.")
    parser.add_argument('--walltime', type = str, default = '0-0:30:0', required = False,\
                        help = "Walltime of each job. Default: 0-0:30:0")
    parser.add_argument('--wait', type = int, default = 1, required = False,\
                        help = "Time interval in seconds between submission of two jobs. Default: 1")
    parser.add_argument('--opts', type = str, required = True,\
                        help = "A string passed to readSimulator.py for arguments")
    
    return parser.parse_args()


def main():
    args = parse_arguments()
    script = os.path.join(args.path, 'readSimulator.py')
    
    if not os.path.exists(script):
        sys.exit("Error: readSimulator.py is not found under the directory " + args.path + ".")
    
    for f in args.input:  # f: absolute or relative path to a FASTA file
        job_id = os.path.splitext(os.path.basename(f))[0]  # remove file name extension
        
        cmd = '#!/bin/bash'
        cmd += '\n#SBATCH -p ' + PARTITION
        cmd += '\n#SBATCH --job-name=' + job_id
        cmd += '\n#SBATCH --ntasks=1'
        cmd += '\n#SBATCH --mem-per-cpu=' + args.mem
        cmd += '\n#SBATCH --time=' + args.walltime
        cmd += '\n\nmodule load ' + MODULE_PYTHON
        cmd += '\nmodule load ' + MODULE_SAMTOOLS
        cmd += '\n\npython %s --input %s %s\n' % (script, f, args.opts)
        
        slurm_filename = job_id + '.slurm'
        with open(slurm_filename, 'w') as slurm_script:
            slurm_script.write(cmd)
            
        os.system('sbatch ' + slurm_filename)  # submit this job
        time.sleep(args.wait)

    
if __name__ == '__main__':
    main()
