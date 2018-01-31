#!/usr/bin/env python

import argparse
import numpy as np
import re

def main():
    
    # Parse arguments.
    parser = argparse.ArgumentParser(description='get contig counts and stats from a multi-fasta file')
    parser.add_argument('-i', '--input', help='indicate input fasta', required=True)
    parser.add_argument('-s', '--show', help='indicate True to print length of each contig', required=False)
    args = parser.parse_args()
    
    # Open FASTA.
    ctgLens = []
    p = re.compile('(>.+\n)|\W')
    N = re.compile('[nN]')
    Ns = 0

    with open(args.input, 'rU') as fastaIn:
    	tempSeq = ''
        ctgid = 'null'
	for line in fastaIn:
    		if line.startswith('>'):
                        if args.show:
                            print(ctgid.strip('>') + '\t' + str(len(tempSeq)))
                            ctgid = line.strip('\n')
    			ctgLens.append(len(tempSeq))
    			Ns += len(N.findall(tempSeq))
                        tempSeq = ''
    		else:
    			tempSeq += p.sub('', line)
    	ctgLens.append(len(tempSeq))
        Ns += len(N.findall(tempSeq))
        if args.show:
            print(ctgid.strip('>') + '\t' + str(len(tempSeq)) + '\n')
    	del(tempSeq)
    #ctgLens.sort()
    ctgLens.pop(0)
    ctgCount = len(ctgLens)
    totalLen = sum(ctgLens)
    ctgMax = max(ctgLens)
    ctgMin = min(ctgLens)
    ctgAvg = int(round(np.mean(ctgLens)))
    ctgMed = int(round(np.median(ctgLens)))
    ctgStd = round(np.std(ctgLens), 2)

    print('for ' + args.input + ':\nnumber of contigs = ' + str(ctgCount) + '\ncombined length = ' + str(totalLen) + 'bp\nmax length = ' + str(ctgMax) + 'bp\nmin length = ' + str(ctgMin) + 'bp\naverage length = ' + str(ctgAvg) + 'bp\nmedian length = ' + str(ctgMed) + 'bp\nSD = ' + str(ctgStd) + '\nNs = ' + str(Ns) + '\n')

if __name__ == '__main__':
    main()
