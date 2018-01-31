#!/usr/bin/env python


import argparse, csv

def main():
    
    # Parse arguments.
    parser = argparse.ArgumentParser(description='convert NLR-parser txt output to bed')
    parser.add_argument('-i', '--input', help='indicate input.txt', required=True)
    parser.add_argument('-o', '--output', help='indicate output.bed', required=True)
    args = parser.parse_args()
    
    file_in = open(args.input, 'rU')
    file_out = open(args.output, 'wb')
    reader_in = csv.reader(file_in, delimiter = '\t')
    count = 0
    
    file_out.write('#track name="NLR-Parser_predicted"\n#itemRgb="Off"\n')
    for row in reader_in:
        if 'SequenceName'in row[0]:
            pass
        else:
            count += 1
            if 'forward' in row[3]:
                strand = '+'
            elif 'reverse' in row[3]:
                strand = '-'
            nuRow = row[0] + '\t' + row[4] + '\t' + row[5] + '\t' + row[0] + '_' + row[2] + '\t0\t' + strand + '\t' + row[4] + '\t' + row[5] + '\t0,0,0\n'
            file_out.write(nuRow)

    file_out.close()
    print('wrote ' + str(count) + ' rows to ' + args.output)


if __name__ == '__main__':
    main()
