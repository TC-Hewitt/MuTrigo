# Copyright (C) 2017 Timothy C. Hewitt - All Rights Reserved
# You may use, distribute and modify this code under the terms of the GNU Public License version 3 (GPLv3)
# You should have recieved a copy of the GPLv3 license with this file. If not, please visit https://github.com/TC-Hewitt/MuTrigo

#!/usr/bin/env python

from __future__ import division
from numpy.random import randint
import argparse, sys, re, csv
csv.field_size_limit(sys.maxsize)

def main():

    # Parse arguments.
    parser = argparse.ArgumentParser(description='SNPlogger will parse an mpileup file and log all SNPs and indels that satisfy parameters. Final tally printed to STDOUT. Outfile is formatted as tab sep fields: <seqid> <position(1based)> <polymorphic-type> <frequency(float)>. Compatible with STDIN.')
    parser.add_argument('-i', '--input', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='indicate input.pileup (leave out if using STDIN). Best if pileups generated with -a/-aa option (samtools > v1.4).')
    parser.add_argument('-o', '--output', help='indicate output file', required=True)
    parser.add_argument('-d', '--mindep', help='set min depth. Only bases with read coverage equal to or above this number are considered for SNP or indel calling (default=10)', default=10, type=int, required=False)
    parser.add_argument('-f', '--minfrq', help='set min frequency of any mismatch at base to call a SNP (default=0.2). Default threshold will call mixed allelic SNVs. Note: Ns in reference not counted for SNPs.', default=0.2, type=float, required=False)
    parser.add_argument('-x', '--idfrq', help='set min frequency of indel to report an indel (default=0.8)', default=0.8, type=float, required=False)
    parser.add_argument('-b', '--blacklist', type=argparse.FileType('r'), help='provide a noisefinder outfile listing contig regions to omit from analysis.', required=False)
    parser.add_argument('-a', '--appendbl', help='indicate a noisefinder outfile to append its contents to SNPlogger out in adjusted format. Or indicate "True" to use same file as in -b/--blacklist (can be useful to include poor coverage/alignment zones in subsequent mutant analysis - <position> field contains start of low coverage or noisy alignment region).', required=False)
    args = parser.parse_args()

    # retrieve contigs from blacklist
    if args.blacklist:
        ctgdict = {}
        listIn = csv.reader(args.blacklist, delimiter = '\t')
        for row in listIn:
            try:
                if row[0] in ctgdict:
                    ctgdict.setdefault(row[0], []).append((int(row[1]),int(row[2])))
                else:
                    ctgdict[row[0]]=[(int(row[1]),int(row[2]))]
            except (IndexError, ValueError):
                continue
        #generates something like: {'contig_1':[(210,510),(1215,3211)],'contig_2':[(123,456),(789,1112),...} 
        print(str(len(ctgdict.keys())) + ' contigs added to blacklist.\n')
    
    # set up counters for depth logging
    above = 0
    below = 0
    Nabove = 0
    Nbelow = 0
    masked = 0

    # set up counters for SNP/del logging
    pileIn = csv.reader(args.input, delimiter = '\t', quoting=csv.QUOTE_NONE)
    fileOut = open(args.output, 'w')
    maxref = 1.0 - args.minfrq
    A = {'T':0,'C':0,'G':0}
    T = {'A':0,'C':0,'G':0}
    C = {'A':0,'T':0,'G':0}
    G = {'A':0,'T':0,'C':0}
    ref = {'A':A, 'T':T, 'C':C, 'G':G}
    bases = {'A':'TCG', 'T':'ACG', 'C':'ATG', 'G':'ATC'}
    indels = 0
    SNPs = 0
    pat1 = re.compile('[atcgn]', re.I)
    pat2 = re.compile('[+-]\d+')

    # open pileup and parse
    current = None
    mask = False
    for row in pileIn:
        if row[0] != current:
            current = row[0]
            if args.blacklist and current in ctgdict:
                mask = True
                zones = ctgdict[current]
            else:
                mask = False
        elif int(row[3]) < args.mindep: # ignore rows if mindep below cutoff
            below += 1
            if row[2] == 'N':
                Nbelow += 1
            continue
        elif not pat1.search(row[4]): # ignore rows if no mismatch present
            above += 1
            if row[2] == 'N':
                Nabove += 1
            continue
        elif len(pat1.findall(row[4]))/int(row[3]) >= args.minfrq:
            if mask == True and any(min <= int(row[1]) <= max for (min,max) in zones): # ignore rows if mask in ON and in zone that is blacklisted.
                masked += 1
                continue
            above += 1
            if row[2] == 'N':
                Nabove += 1
            mmatches = ''.join(pat1.findall(row[4])).upper()
            truPos = False
            dep = int(row[3])
            if not pat2.search(row[4]):
                for k in bases.keys():
                    if row[2] == k:
                        for b in bases[k]:
                            freq = mmatches.count(b)/dep
                            if freq >= args.minfrq:
                                if truPos == False:
                                    SNPs += 1
                                truPos = True
                                ref[k][b] += 1
                                fileOut.write(row[0] + '\t' + row[1] + '\t' + k + '>' + b + '\t' + str(round(freq,3)) + '\n')
                            else:
                                continue
                            if freq > maxref:
                                break
                        break
            else:
                InDel = pat2.findall(row[4])
                freq = len(InDel)/dep
                if freq >= args.idfrq:
                    fileOut.write(row[0] + '\t' + row[1] + '\tindel>' + ','.join(list(set(InDel))) + '\t' + str(round(freq,3)) + '\n')
                    indels += 1
        else: # ignore rows if overall mismatch rate below cutoff
            above += 1
            if row[2] == 'N':
                Nabove += 1
            continue

    # append contents of noisefinder if indicated by -b
    if args.appendbl:
        lowcov = 0
        noisy = 0
        if args.appendbl in ['T', 't', 'True', 'true', 'TRUE']:
            args.blacklist.seek(0)
            for row in listIn:
                randVal = randint(0, 100, 1) # a random number is added to the start coord so that SNPtracker won't disregard noisy/lowcov features with identical starts in multiple mutants
                try:
                    randStart = int(row[1]) + randVal[0]
                    if 'xxx' in row[4]:
                        lowcov += 1
                        fileOut.write(row[0] + '\t' + str(randStart) + '\tlowcov\tNaN\n')
                    elif float(row[4]):
                        noisy += 1
                        fileOut.write(row[0] + '\t' + str(randStart) + '\tnoisy\tNaN\n')
                except (IndexError, ValueError):
                    continue
        else:
            appendIn = open(args.appendbl, 'rU')
            appendOut = csv.reader(appendIn, delimiter = '\t')
            for row in appendOut:
                randVal = randint(0, 100, 1)
                try:
                    randStart = int(row[1]) + randVal[0]
                    if 'xxx' in row[4]:
                        lowcov += 1
                        fileOut.write(row[0] + '\t' + str(randStart) + '\tlowcov\tNaN\n')
                    elif float(row[4]):
                        noisy += 1
                        fileOut.write(row[0] + '\t' + str(randStart) + '\tnoisy\tNaN\n')
                except (IndexError, ValueError):
                    continue
    fileOut.close()

    total = above + below + masked
    print('SNP positions detected: ' + str(SNPs) + '\n<type>\t<occurences>\n')
    for i in A:
        print('A>' + i + ':\t' + str(A[i]))
    print('')
    for i in T:
        print('T>' + i + ':\t' + str(T[i]))
    print('')
    for i in C:
        print('C>' + i + ':\t' + str(C[i]))
    print('')
    for i in G:
        print('G>' + i + ':\t' + str(G[i]))
    print('\nindels' + ':\t' + str(indels) + '\n')
    if args.appendbl:
        print('appended ' + str(lowcov) + ' low coverage regions and '+ str(noisy) + ' noisy alignment regions to output.\n')
    print(str(args.input) + '\ndepth cutoff: ' + str(args.mindep) + '\nbp total=' + str(total) + '\nbp above=' + str(above) + ' (' + str(Nabove) + ' Ns)' + '\nbp below=' + str(below) + ' (' + str(Nbelow) + ' Ns)' + '\nSNPs masked=' + str(masked) + '\n')

if __name__ == '__main__':
    main()
