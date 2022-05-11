# Copyright (C) 2017 Timothy C. Hewitt - All Rights Reserved
# You may use, distribute and modify this code under the terms of the GNU Public License version 3 (GPLv3)
# You should have recieved a copy of the GPLv3 license with this file. If not, please visit https://github.com/TC-Hewitt/MuTrigo

#!/usr/bin/env python

# takes arbitrary number of arguments for either WT or mutant logs
# WTs are read into memory as sets (retaining seqID and coordinate) which are then unioned to a single set
# for each mutant log, it is read into memory as a set (retaining seqID and coordinate) and then has the WT set subtracted from it
# each of the reduced mutant lists then has the coordinate suffix stripped from it so only the seqID remains
# each combination of all possible subset of mutant sets are generated
# for each combination subset of mutant sets, the intersection is found and written to a summary report as well as a detailed report if option given

import argparse, sys, csv, math, re, itertools
csv.field_size_limit(sys.maxsize)

def main():

    # Parse arguments.
    parser = argparse.ArgumentParser(description='finds sequence IDs/regions with coinciding polymorphic features across multiple SNPlogger generated files')
    parser.add_argument('-w', '--wildtype', help='indicate space sep list of logfiles whose features to mask from mutant logfiles', nargs='*', required=False)
    parser.add_argument('-m', '--mutant', help='indicate space sep list of mutant logfiles', nargs='*', required=True)
    parser.add_argument('-o', '--output', help='indicate prefix only of output html(s)', required=False, default='SNPtracker')
    parser.add_argument('-s', '--select', help='selective by polymorphism type. Indicate space sep list of types to only include in analysis (default includes all). Accepted strings are any base change in the form N\>N (eg. "C\>T"), "indel", "lowcov", "noisy", "any" (any N\>N)', nargs='*', type=str, required=False)
    parser.add_argument('-f', '--filter', help='filter by SNV frequency. Indicate min frequency to include in analysis (default=0.8). Entries with "NaN" included by default', type=float, required=False, default=0.8)
    parser.add_argument('-v', '--verbose', help='indicate True to also generate detailed reports (incl. polymorphic type and coordinate) for each subset number of mutants in addition to default summary html', type=str, required=False)
    parser.add_argument('-p', '--proximal', help='indicate window size. Instead of finding features that coincide on a particular contig, SNPtracker will find features that reside close to each other within a user defined window size (min=1000 bases). Suitable for large scaffolds or pseudomolecules', type=int, required=False)
    parser.add_argument('-n', '--min', help='set min number of mutants to consider. Otherwise all mutant subsets >=2 are analysed', type=int, required=False, default=2)
    parser.add_argument('-t', '--tolerate', help='set max number of mutants to tolerate with polymorphisms in identical positions for a given discovery (default none)', type=int, required=False, default=1)
    args = parser.parse_args()

    if len(args.mutant) < 2:
        sys.exit("--mutant needs at least 2 arguments!")

    mVarD = {} # mutant var names paired with original filename
    mRedD = {} # mutant var names paired with redundant set of their concat features and position
    mSetD = {} # mutant var names paired with redundancy removed set of their concat features and position
    mRawD = {} # mutant var names paired with nested dict that pairs contig id with list of tuples containing all features
    mSubD = {} # subsets (N>=2) paired with the combos of mutants sub(N)
    mMask = {} # mutant var names paired with set of concat features and positions to mask based on args.filter|select

    if args.select:
        if 'any' in args.select:
            args.select.remove('any')
            args.select.append('\\w>\\w')
        features = re.compile('(' + ')|('.join(args.select) + ')')

    if args.verbose in ['T', 't', 'True', 'true', 'TRUE']:
        verbose = True
    else:
        verbose = False

    if args.proximal:
        if args.proximal < 1000: # window size set to its allowed minumum if below that
            print('--proximal window size too small. Setting to 1000!')
            args.proximal = 1000
        tile = args.proximal//4 # sets an arbitrary tile size of 1/4 the given window size to increment by

    if args.tolerate < 1 or args.tolerate > len(args.mutant):
        print('cannot use number given for --tolerate. Setting to default!')
        args.tolerate = 1

    wSet = set([])
    if args.wildtype:
        if len(args.wildtype) > 1:
            for i in range(len(args.wildtype)):
                temp = open(args.wildtype[i], 'r')
                wLog = csv.reader(temp, delimiter = '\t', quoting=csv.QUOTE_NONE)
                for row in wLog:
                    wSet.add(row[0] + ' ' + row[1])
            temp.close()
        elif len(args.wildtype) == 1:
            temp = open(args.wildtype[0], 'r')
            wLog = csv.reader(temp, delimiter = '\t', quoting=csv.QUOTE_NONE)
            for row in wLog:
                wSet.add(row[0] + ' ' + row[1])
            temp.close()

    for i in range(len(args.mutant)): # iterate over mutant files
        mVar = 'm'+str(i)
        mVarD[mVar] = args.mutant[i]
        mSet = set([])
        mMask[mVar] = set([]) # used only if verbose combined with filter or select
        with open(args.mutant[i], 'r') as temp:
            mLog = csv.reader(temp, delimiter = '\t', quoting=csv.QUOTE_NONE)
            if args.select and args.filter:
                for row in mLog:
                    if re.match(features, row[2]) and (float(row[3]) >= args.filter or math.isnan(float(row[3]))):
                        mSet.add(row[0] + ' ' + row[1])
                    elif verbose == True or args.proximal:
                        mMask[mVar].add(row[0] + ' ' + row[1])
            elif args.select or args.filter:
                for row in mLog:
                    keep = False
                    if args.filter and (float(row[3]) >= args.filter or math.isnan(float(row[3]))):
                        keep = True
                    elif args.select and re.match(features, row[2]):
                        keep = True
                    if keep == True:
                        mSet.add(row[0] + ' ' + row[1])
                    elif verbose == True or args.proximal:
                        mMask[mVar].add(row[0] + ' ' + row[1])
            else:
                for row in mLog:
                    mSet.add(row[0] + ' ' + row[1])
            if args.wildtype:
                rSet = mSet - wSet
                print(str(len(mSet) - len(rSet)) + ' features shared with wildype(s) masked from ' + mVarD[mVar] + ' after selection/filtering.')
                del(mSet)
                mRedD[mVar] = rSet
                del(rSet)
            else:
                mRedD[mVar] = mSet
                del(mSet)

    # iterate over subset combos from N(mutants) to --tolerate to remove identical SNPs
    mNums = []
    for n in range(args.tolerate+1,len(mVarD)+1):
        mSub = 'N'+str(n)
        mNums.append(mSub)
        mSubD[mSub] = []
        for subset in itertools.combinations(mVarD.keys(), n):
            mSubD[mSub].append(subset)
    mNums = mNums[::-1]
    rGlobal = set([]) # set of globally redundant positions found during runtime
    for n in mNums: # iterate over n mutants
        if args.tolerate == len(args.mutant):
            break
        for s in mSubD[n]: # iterate over combinations of n mutants
            rTemp = []
            for var in s:
                rTemp.append(mRedD[var])
            rInter = set.intersection(*rTemp) # intersection of mutant subset (identical/redundant positions across n mutants)
            rGlobal.update(rInter)
    for mVar in mVarD:
        rmdup = mRedD[mVar] - rGlobal # remove any globally redundant positions from each mutant
        mSetD[mVar] = set([re.sub('\s\d+$','',x) for x in list(rmdup)]) # removes space and coordinate leaving only seq IDs in set
        if verbose == True or args.proximal:
            mRawD[mVar] = dict.fromkeys(mSetD[mVar], None)
            for k in mRawD[mVar].keys():
                mRawD[mVar][k]=[]
            temp = open(mVarD[mVar], 'r')
            mLog = csv.reader(temp, delimiter = '\t', quoting=csv.QUOTE_NONE)
            maskTotal = wSet | mMask[mVar] | rGlobal
            for row in mLog:
                if row[0] in mRawD[mVar] and ' '.join([row[0],row[1]]) not in maskTotal:
                    mRawD[mVar][row[0]].append(tuple(row[1:]))
            temp.close()

    del(wSet)
    del(mMask)
    del(mRedD)
    
    nameOut = str(args.output) + '_summary.html'
    summaryOut = open(nameOut, 'w+')
    summaryOut.write('<!DOCTYPE html>\n<html>\n<h1>summary</h1>\n<h3>parameters</h3>\n')
    if args.wildtype:
        summaryOut.write('<p>\nwildtypes: ' + ', '.join(args.wildtype))
    else:
        summaryOut.write('<p>\nwildtypes: NA')
    summaryOut.write('<br>\nmutants: ' + ', '.join(args.mutant))
    if args.select:
        if '\\w>\\w' in args.select:
            args.select.remove('\\w>\\w')
            args.select.append('any')
        summaryOut.write('<br>\nselected: ' + ', '.join(args.select))
    else:
        summaryOut.write('<br>\nselected: NA')
    if args.filter:
        summaryOut.write('<br>\nfiltered: ' + str(args.filter))
    else:
        summaryOut.write('<br>\nfiltered: NA')
    if args.proximal:
        summaryOut.write('<br>\nproximal: ON, window: ' + str(args.proximal))
    else:
        summaryOut.write('<br>\nproximal: OFF')
    if args.tolerate != 1:
        summaryOut.write('<br>\ntolerate: ' + str(args.tolerate) + '\n</p>\n')
    else:
        summaryOut.write('<br>\ntolerate: none\n</p>\n')

    del(mNums[:])
    mSubD.clear()
    for n in range(args.min,len(mSetD)+1):
        mSub = 'N'+str(n)
        mNums.append(mSub)
        mSubD[mSub] = []
        for subset in itertools.combinations(mSetD.keys(), n):
            mSubD[mSub].append(subset)
    mNums = mNums[::-1]
    sGlobal = set([]) # set of all contigs already found during runtime
    for n in mNums: # iterate over n mutants
        nHits = 0
        nInt = int(n.strip('N'))
        summaryOut.write('\n<h3>polymorphic in ' + str(nInt) + ' mutants</h3>\n<p>\n')
        if verbose == True: # open report file to write to if verbose true
            vnameOut = str(args.output) + '_' + n + '_report.html'
            verboseOut = open(vnameOut, 'w+')
            verboseOut.write('<!DOCTYPE html>\n<html>\n<body>\n<h1>polymorphic in ' + str(nInt) + ' mutants</h1>\n')
        for s in mSubD[n]: # iterate over combinations of n mutants
            sTemp = []
            mNames = [mVarD[var] for var in s]
            for var in s:
                sTemp.append(mSetD[var])
            sInter = set.intersection(*sTemp) - sGlobal # diff of intersection of mutant subset minus contigs already found (prevent duplication if promixmal off)
            if not args.proximal:
                sGlobal.update(sInter)
            if args.proximal and len(sInter) != 0:
                for seq in sInter: # within contig testing
                    coords = {}
                    allVals = []
                    for var in s: # retrieve coords from mRawD[var], setup sliding window loop testing each time below
                        try:
                            coords[var] = [int(feature[0]) for feature in mRawD[var][seq]]
                        except BaseException as err:
                            print(err.message)
                            continue
                    for coord in coords.values():
                        allVals = allVals + coord
                    maxVal = max(allVals)
                    lowLim = min(allVals)
                    if (maxVal - lowLim) <= args.proximal: # write out all mutant features if min to max range is lower than given proximal range
                        nHits += 1
                        summaryOut.write(seq + ':' + str(lowLim) + '-' + str(maxVal) + ' <----- (' + ', '.join(mNames) + ')<br>\n')
                        if verbose == True:
                            verboseOut.write('<h3>' + seq + ':' + str(lowLim) + '-' + str(maxVal) + '</h3>\n<p>\n')
                            for var in s:
                                if seq in mRawD[var]:
                                    verboseOut.write(mVarD[var] + ' ' + str(mRawD[var][seq]).replace('\'','') + '<br>\n')
                            verboseOut.write('</p>\n')
                        continue
                    uppLim = lowLim + args.proximal
                    zones = set([])
                    while lowLim < (maxVal - tile): # within window testing of coords for each mutant: if at least nInt independent features (tally == nInt) given that any(coords per var fall within window), post the coords of lower feature and upper feature in that window to zones (zone set object in case duplicates due to tiling)
                        tally = 0
                        inWindow = []
                        for var in coords:
                            varHits = [coord for coord in coords[var] if lowLim <= coord <= uppLim]
                            if len(varHits) > 0:
                                tally += 1
                                inWindow = inWindow + varHits
                        if tally == nInt:
                            zmin = min(inWindow)
                            zmax = max(inWindow)
                            zones.add((zmin,zmax))
                        lowLim += tile
                        uppLim += tile
                    seqHits = []
                    if len(zones) != 0:
                        nHits += len(zones)
                        for zone in zones:
                            record = seq + ':' + str(zone[0]) + '-' + str(zone[1])
                            seqHits.append(record) # in form eg. "contig_888:1500-3500, contig_888:7000-9000, contig_901:1-2000 (mut1.log, mut3.log, mut5.log)"
                            if verbose == True:
                                verboseOut.write('<h3>' + record + '</h3>\n<p>\n')
                                for var in s:
                                    try:
                                        for feature in mRawD[var][seq]:
                                            if zone[0] <= int(feature[0]) <= zone[1]:
                                                verboseOut.write(mVarD[var] + ' [' + str(feature).replace('\'','') + ']<br>\n')
                                    except BaseException as err:
                                        print(err.message)
                                        continue
                                verboseOut.write('</p>\n')
                        summaryOut.write(', '.join(seqHits) + ' <----- (' + ', '.join(mNames) + ')<br>\n')                
            elif len(sInter) != 0:
                nHits += len(sInter)
                summaryOut.write(', '.join(sInter) + ' <----- (' + ', '.join(mNames) + ')<br>\n')
                if verbose == True: # write to corresponding verbose file
                    for seq in sInter:
                        verboseOut.write('<h3>' + seq + '</h3>\n<p>\n')
                        for var in s:
                            if seq in mRawD[var]:
                                verboseOut.write(mVarD[var] + ' ' + str(mRawD[var][seq]).replace('\'','') + '<br>\n')
                        verboseOut.write('</p>\n')
            else:
                continue
        summaryOut.write('</p>\n')
        print('found across ' + str(nInt) + ' mutants: ' + str(nHits))
        if verbose == True:
            verboseOut.write('</body>\n</html>\n')
            verboseOut.close()
    summaryOut.write('</body>\n</html>\n')
    print('done.')

if __name__ == '__main__':
    main()
