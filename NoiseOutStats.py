#!/usr/bin/env python

import argparse, csv
import numpy as np

def main():
    
    # Parse arguments.
    parser = argparse.ArgumentParser(description='get counts of noisy regions and lowcov regions from Noisefinder out file')
    parser.add_argument('-i', '--input', help='indicate input file', required=True)
    parser.add_argument('-r', '--reflen', help='indicate total length of reference fasta used for mapping', type=float, required=False)
    args = parser.parse_args()
    
    # establish params
    NoiseLens = []
    lowcovLens = []
    binLims = (0, 500, 1000, 2000, 5000, 8000, 12000, 20000, 40000, 60000, 100000, float('inf'))
    binCountsNoise = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
    binCountsLowcov = {0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0}
    binDex = ('0-0.5Kb', '0.5-1Kb', '1-2Kb', '2-5Kb', '5-8Kb', '8-12Kb', '12-20Kb', '20-40Kb', '40-60Kb', '60-100Kb', '100Kb+')
    
    def isfloat(value):
        try:
            float(value)
            return True
        except ValueError:
            return False

    # parse input
    Input = open(args.input, 'r')
    Intab = csv.reader(Input, delimiter = '\t', quoting=csv.QUOTE_NONE)
    for row in Intab:
        try:
            c3Value = int(row[3])
            if isfloat(row[4]) == True:
                i = 0
                while i < 11:
                    if binLims[i] < c3Value <= binLims[i+1]:
                        binCountsNoise[i] += 1
                        break
                    i += 1
                NoiseLens.append(c3Value)
            elif 'xxx' in row[4]:
                i = 0
                while i < 11:
                    if binLims[i] < c3Value <= binLims[i+1]:
                        binCountsLowcov[i] += 1
                        break
                    i += 1
                lowcovLens.append(c3Value)
        except (IndexError, ValueError):
            continue
    
    NoiseLens.sort(reverse=False)
    lowcovLens.sort(reverse=False)

    # get noise stats
    if len(NoiseLens) > 0:
        NoiseCount = len(NoiseLens)
        NoiseLenTotal = sum(NoiseLens)
        if args.reflen:
            NoisePcntRef = (NoiseLenTotal*100)/args.reflen
        else:
            NoisePcntRef = "NaN"
        NoiseLenMax = max(NoiseLens)
        NoiseLenMin = min(NoiseLens)
        NoiseLenMed = int(round(np.median(NoiseLens)))
        NoiseLenQ1 = int(round(np.percentile(NoiseLens, 25)))
        NoiseLenQ3 = int(round(np.percentile(NoiseLens, 75)))
        NoiseLenMean = int(round(np.mean(NoiseLens)))
        NoiseLenStd = round(np.std(NoiseLens), 2)
    else:
        NoiseCount = "NaN"
        NoiseLenTotal = "NaN"
        NoisePcntRef = "NaN"
        NoiseLenMax = "NaN"
        NoiseLenMin = "NaN"
        NoiseLenMed = "NaN"
        NoiseLenQ1 = "NaN"
        NoiseLenQ3 = "NaN"
        NoiseLenMean = "NaN"
        NoiseLenStd = "NaN"

    # get lowcov stats
    if len(lowcovLens) > 0:
        lowcovCount = len(lowcovLens)
        lowcovLenTotal = sum(lowcovLens)
        if args.reflen:
            lowcovPcntRef = (lowcovLenTotal*100)/args.reflen
        else:
            lowcovPcntRef = "NaN"
        lowcovLenMax = max(lowcovLens)
        lowcovLenMin = min(lowcovLens)
        lowcovLenMed = int(round(np.median(lowcovLens)))
        lowcovLenQ1 = int(round(np.percentile(lowcovLens, 25)))
        lowcovLenQ3 = int(round(np.percentile(lowcovLens, 75)))
        lowcovLenMean = int(round(np.mean(lowcovLens)))
        lowcovLenStd = round(np.std(lowcovLens), 2)
    else:
        lowcovCount = "NaN"
        lowcovLenTotal = "NaN"
        lowcovPcntRef = "NaN"
        lowcovLenMax = "NaN"
        lowcovLenMin = "NaN"
        lowcovLenMed = "NaN"
        lowcovLenQ1 = "NaN"
        lowcovLenQ3 = "NaN"
        lowcovLenMean = "NaN"
        lowcovLenStd = "NaN"

    # output table
    print('#<lengths>\t<noisy regions>\t<lowcov regions>')
    print('#Number\t' + str(NoiseCount) + '\t' + str(lowcovCount) +\
    '\n#Total\t' + str(NoiseLenTotal) + '\t' + str(lowcovLenTotal) +\
    '\n#Reference%\t' + str(NoisePcntRef) + '\t' + str(lowcovPcntRef) +\
    '\n#Max\t' + str(NoiseLenMax) + '\t' + str(lowcovLenMax) +\
    '\n#Min\t' + str(NoiseLenMin) + '\t' + str(lowcovLenMin) +\
    '\n#Median\t' + str(NoiseLenMed) + '\t' + str(lowcovLenMed) +\
    '\n#Q1\t' + str(NoiseLenQ1) + '\t' + str(lowcovLenQ1) +\
    '\n#Q3\t' + str(NoiseLenQ3) + '\t' + str(lowcovLenQ3) +\
    '\n#Mean\t' + str(NoiseLenMean) + '\t' + str(lowcovLenMean) +\
    '\n#SD\t' + str(NoiseLenStd) + '\t' + str(lowcovLenStd))
    for i in range(11):
        print('#' + str(binDex[i]) + '\t' + str(binCountsNoise[i]) + '\t' + str(binCountsLowcov[i]))
        
if __name__ == '__main__':
    main()

