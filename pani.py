#!/bin/env python3
#./parallel_ani.py -o output.out -t 4 genome1.fasta genome2.fasta genome3.fasta genome4.fasta genome5.fasta genome6.fasta genome7.fasta genome8.fasta genome9.fasta genome10.fasta 
#Ran command above in ~3:30 seconds, 4 threads with all 10 files. Produced output.out file with matrix
from multiprocessing import Pool
import itertools
import os
import argparse
import subprocess
from pathlib import Path
import re

def cli_parse():
    parser = argparse.ArgumentParser(description='Parallel ANI calculation')
    parser.add_argument('-o',
                        help='This is output file',
                        default='output.txt')
    parser.add_argument('-t',
                        help='This is number of threads',
                        type=int,
                        default=1,
                        choices=range(1,5))
    parser.add_argument('strings',
                        metavar='N',
                        nargs='+',
                        help='Microbial FASTA files')
    return parser.parse_args()

def cleanup():
    delete = {'.mcoords','.mdelta','.qdiff','.snps'
		,'.rdiff','.report','.1delta','.delta','.1coords'}
    for p in Path.cwd().iterdir():
        if p.suffix in delete:
            p.unlink()

def data(args_arg):
    return list(itertools.combinations(args_arg, r=2))

def parallel(num):
    print(f'Process {os.getpid()} working record {num}')
    temp = str(os.getpid())
    temp1 = temp + '.report'
    v = subprocess.check_output(['dnadiff',
                                 '-p',
                                 temp,
                                 num[0],
                                 num[1]])
    with open(temp1, mode='r') as file:
        for line in file:
            if line.startswith('AvgIdentity'):
                var = (frozenset({num[0],num[1]}), line.split()[1])
                break
    return var 

def file_write(args,identity):
    with open(args.o, mode='w') as file:
        for row, row_val in enumerate(args.strings):
            for col, col_val in enumerate(args.strings):
                if row == col:
                    file.write('100\t')
                elif col == (len(args.strings) - 1):
                    text = str(identity[frozenset({row_val,col_val})]) + '\n'
                    file.write(text)
                else:
                    text = str(identity[frozenset({row_val,col_val})]) + '\t'
                    file.write(text)

def main():
    args = cli_parse()
    formatted = data(args.strings)
    pool = Pool(args.t)
    identities = dict(pool.map(parallel,formatted))
    pool.close()
    pool.join()
    file_write(args,identities)
    cleanup()

if __name__ == "__main__":
    main()

