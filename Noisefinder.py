# Copyright (C) 2017 Timothy C. Hewitt - All Rights Reserved
# You may use, distribute and modify this code under the terms of the GNU Public License version 3 (GPLv3)
# You should have recieved a copy of the GPLv3 license with this file. If not, please write to timcharleshewitt@gmail.com or visit https://github.com/TC-Hewitt/MuTrigo

#!/usr/bin/env python

from __future__ import division
import argparse, csv, sys, re
csv.field_size_limit(sys.maxsize)

def main():

    # Parse arguments.
    parser = argparse.ArgumentParser(description='Regions rich in mismatches/poor coverage after read alignment can often signify misalignment or mixed alignment due to allelism, polyploidy, or large deletions. Given a pileup file, noisefinder reports regions containing a density of SNVs above a user defined threshold over a given min length and min read depth (prints to STDOUT).')
    parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='indicate input pileup. Leave out option if piping from STDIN. Best if pileups generated with -a/-aa option (samtools > v1.4)')
    parser.add_argument('-d', '--mindep', help='set min depth. Only bases with read coverage equal to or above this number are considered for SNV calling (default=5)', default=5, type=int, required=False)
    parser.add_argument('-l', '--minlen', help='set min length. Only compute SNV frequency for regions above this length (default=300)', default=300, type=int, required=False)
    parser.add_argument('-c', '--regnf', help='set min density (frequency over region) of SNVs. Only report regions, contigs having a SNV density higher than or equal to this (default=0.005 aka 1/200 bases)', default=0.005, type=float, required=False)
    parser.add_argument('-b', '--basef', help='set min frequency of mismatch at base to call a SNV (default=0.2)', default=0.2, type=float, required=False)
    parser.add_argument('-a', '--addlc', help='indicate min length of regions below depth cutoff to include in final output (these are not SNV counted but marked with "xxx" in last field)', type=int, required=False)
    args = parser.parse_args()

    # open pileup
    pileIn = csv.reader(args.infile, delimiter = '\t', quoting=csv.QUOTE_NONE)
    print('#parsing pileup...\n#\n#<seq_id>\t<start>\t<end>\t<length>\t<SNV_density>')
    # global metrics
    scontigs = 0
    regions = 0
    # regional metrics
    snvs = 0
    bplen = 0
    contig = ''
    start = '' 
    end = ''
    rpg = 0 # stands for "regions per contig" that satisfy criteria 
    # set switch, parse rows (switch is on "True" when a candidate region is being counted)
    switch = False
    incLDR = False # include low depth regions. Set to true when -a given an integer argument
    if args.addlc:
        if args.addlc < 200:
            sys.exit("option --addlc does not accept lengths less than 200.")
        incLDR = True
        lenLDR = 0
        sttLDR = 0
        endLDR = ''
        ctgLDR = ''
    pat = re.compile('[atcgn]', re.I)
    for row in pileIn:
        if int(row[3]) < args.mindep and switch == False: # ignore rows if mindep below cutoff and switch is off
            if incLDR:
                if row[0] == ctgLDR: # count rows if in low depth region and --addlc option is on
                    lenLDR += 1
                else: # print current LDR stats if seqID (row 0) changes and reinitialise LDR stats for new seqID
                    endLDR = int(sttLDR) + lenLDR
                    if lenLDR >= args.addlc:
                        print(ctgLDR + '\t' + sttLDR + '\t' + str(endLDR) + '\t' + str(lenLDR) + '\txxx')
                    ctgLDR = row[0]
                    sttLDR = row[1]
                    lenLDR = 1
            continue
        elif int(row[3]) >= args.mindep and switch == False: # initialise new candidate region if mindep rises above cutoff and switch previously off
            bplen = 1
            switch = True
            start = row[1]
            snvs = 0
            if row[0] != contig:
                contig = row[0]
                rpg = 0
            if len(pat.findall(row[4]))/int(row[3]) >= args.basef:
                snvs += 1
            if incLDR: # print current LDR stats if --addlc is on
                endLDR = int(sttLDR) + lenLDR
                if lenLDR >= args.addlc:
                    print(ctgLDR + '\t' + sttLDR + '\t' + str(endLDR) + '\t' + str(lenLDR) + '\txxx')
        elif int(row[3]) >= args.mindep and row[0] == contig and switch == True: # while switch is on and mindep stays above cutoff, rows will be SNV tested
            bplen += 1
            if len(pat.findall(row[4]))/int(row[3]) >= args.basef:
                snvs += 1
        elif int(row[3]) < args.mindep and row[0] == contig and switch == True: # if mindep drops below cutoff while switch is on, signals end of region and prints results if all criteria satisfied. Metrics reset
            if bplen > args.minlen and snvs/bplen >= args.regnf:
                regions += 1
                end = row[1]
                rpg += 1
                if rpg == 1:
                    scontigs += 1
                print(contig + '\t' + start + '\t' + end + '\t' + str(bplen) + '\t' + str(round(snvs/bplen, 3)))
            bplen = 0
            snvs = 0
            switch = False
            if incLDR: # reinitialise LDR stats for new region if --addlc is on
                ctgLDR = row[0]
                sttLDR = row[1]
                lenLDR = 1
        elif str(row[0]) != contig and switch == True: # if contig id changes while switch is on, signals end of region and prints results if all criteria satisfied. Metrics reset
            if bplen > args.minlen and snvs/bplen >= args.regnf:
                end = int(start) + bplen
                regions += 1
                rpg += 1
                if rpg == 1:
                    scontigs += 1
                rpg = 0
                print(contig + '\t' + start + '\t' + str(end) + '\t' + str(bplen) + '\t' + str(round(snvs/bplen, 3)))
                contig = row[0]
            if int(row[3]) >= args.mindep: # if mindep of first base in new contig satisfies cutoff, new candidate region initialised
                start = row[1]
                bplen = 1
                if len(pat.findall(row[4]))/int(row[3]) >= args.basef:
                    snvs = 1
            else: # metrics reset if mindep below cutoff
                bplen = 0
                snvs = 0
                switch = False
                if incLDR: # reinitialise LDR stats for new region if --addlc is on
                    ctgLDR = row[0]
                    sttLDR = row[1]
                    lenLDR = 1
    if switch == True: # prints final region results if switch still on when end of file reached and all criteria satisfied
        if bplen > args.minlen and snvs/bplen >= args.regnf:
            regions += 1
            end = int(start) + bplen
            rpg += 1
            if rpg == 1:
                scontigs += 1
            print(contig + '\t' + start + '\t' + str(end) + '\t' + str(bplen) + '\t' + str(round(snvs/bplen, 3)))
    elif incLDR: # prints final LDR results if switch off when end of file reached and --addlc is on
        endLDR = int(sttLDR) + lenLDR
        if lenLDR >= args.addlc:
            print(ctgLDR + '\t' + sttLDR + '\t' + str(endLDR) + '\t' + str(lenLDR) + '\txxx')
    print('########\n#in ' + str(scontigs) + ' contigs, found ' + str(regions) + ' regions of length ' + str(args.minlen) + ' or more containing a SNV density of at least ' + str(args.regnf) + ' with a min frequency of ' + str(args.basef) + ' to call as SNV.')

if __name__ == '__main__':
    main()
