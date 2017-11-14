tools and pipelines that help with deep sequencing analysis

# Variant discovery
this pipeline was optimised for mutagen induced SNP discovery in sequence data from exome capture experiments. However, flexible options allow wide application for variant analysis in any deep sequencing data.

written for Python 2.7 - certain features may not be compatible with Python 3

<Noisefinder.pyc>

    Description: 
        Regions rich in mismatches/poor coverage after read alignment can often signify misalignment or mixed alignment due to allelism, polyploidy, or presence/absemce variation. Given a pileup file, noisefinder reports regions containing a density of SNVs above a user defined threshold over a given min length and min read depth (prints to STDOUT).
    Options:
        '-i', '--infile', 'indicate input pileup. Leave out option if piping from STDIN. Best if pileups generated with -a/-aa option (samtools > v1.4)'
        '-d', '--mindep', 'set min depth. Only bases with read coverage equal to or above this number are considered for SNV calling (default=5)'
        '-l', '--minlen', 'set min length. Only compute SNV frequency for regions above this length (default=300)'
        '-c', '--regnf', 'set min density (frequency over region) of SNVs. Only report regions, contigs having a SNV density higher than or equal to this (default=0.005 aka 1/200 bases)'
        '-b', '--basef', 'set min frequency of mismatch at base to call a SNV (default=0.2)'
        '-a', '--addlc', 'indicate True to include regions below depth cutoff in final output (these are not SNV counted but marked with "xxx" in last field)'
    Library dependencies:
        '__future__', 'argparse', 'csv', 'sys', 're'
    
<SNPlogger.pyc>

    Description:
        SNPlogger will parse an mpileup file and log all SNPs and indels that satisfy parameters. Final tally printed to STDOUT. Outfile is formatted as tab sep fields: <seqid> <position(1based)> <polymorphic-type> <frequency(float)>. Compatible with STDIN.
    Options:
        '-i', '--input', 'indicate input.pileup (leave out if using STDIN). Best if pileups generated with -a/-aa option (samtools > v1.4).')
        '-o', '--output', 'indicate output file'
        '-d', '--mindep', 'set min depth. Only bases with read coverage equal to or above this number are considered for SNP or indel calling (default=5)'
        '-f', '--minfrq', 'set min frequency of any mismatch at base to call a SNP (default=0.2). Default threshold will call mixed allelic SNVs. Note: Ns in reference not counted for SNPs.'
        '-x', '--idfrq', 'set min frequency of indel to report an indel (default=0.8)'
        '-b', '--blacklist', 'provide a noisefinder outfile listing contig regions to omit from analysis.'
        '-a', '--appendbl', 'indicate a noisefinder outfile to append its contents to SNPlogger out in adjusted format. Or indicate "True" to use same file as in -b/--blacklist (can be useful to include poor coverage/alignment zones in subsequent mutant analysis - <position> field contains start of low coverage or noisy alignment region rounded to nearest hundred).'
    Library dependencies:
        '__future__', 'argparse', 'csv', 'sys', 're'
    
<SNPtracker.pyc>

    Description:
        Finds sequence IDs/regions with coinciding variant features across multiple SNPlogger generated files.
    Options:
        '-w', '--wildtype', 'indicate space sep list of logfiles whose features to mask from mutant logfiles', nargs='*', required=False)
        '-m', '--mutant', 'indicate space sep list of mutant logfiles'
        '-o', '--output', 'indicate prefix name of output report(s) (default="SNPtracker")'
        '-s', '--select', 'selective by polymorphism type. Indicate space sep list of types to only include in analysis (default includes all). Accepted strings are any base change in the form N\>N (eg. "C\>T"), "indel", "lowcov", "noisy", "any" (any N\>N)'
        '-f', '--filter', 'filter by SNV frequency. Indicate min frequency to include in analysis (default=0.8). Entries with "NaN" included by default'
        '-v', '--verbose', 'indicate True to also generate detailed reports (incl. polymorphic type and coordinate) for each subset number of mutants in addition to default summary report'
        '-p', '--proximal', 'indicate window size. Instead of finding features that coincide on a particular contig, SNPtracker will find features that reside close to each other within a user defined window size (min=1000 bases). Suitable for assemblies with large scaffolds'
        '-n', '--min', 'set min number of mutants to consider. Otherwise all mutant subsets >=2 are analysed'
        '-t', '--tolerate', 'set max number of mutants to tolerate with polymorphisms in identical positions for a given discovery (default none)'
    Library dependencies:
        'argparse', 'sys', 'csv', 'math', 're', 'itertools'
        
Example Workflow

    coming soon...

# Miscellaneous

<blast_filterV2.pyc>
    
    Description:
        filter and sort a blast output (outfmt 6 or 7) based on the cutoffs you want for each parameter.
    Options:
        '-i', '--input', 'indicate input file (as blast tab outfmt 6 or 7)'
        '-o', '--output', 'indicate output file'
        '-p', '--pcntid', 'set min percent id', default=0, 
        '-a', '--alnlen', 'set min alignment length', default=0
        '-m', '--msmtch', 'set max mismatches', default=infinite
        '-g', '--gpopen', 'set max gap opens', default=infinite
        '-qs', '--qstart', 'set min query seq start', default=0
        '-qe', '--q_end', 'set max query seq end', default=infinite
        '-ss', '--sstart', 'set min subject seq start', default=0
        '-se', '--s_end', 'set max subject seq end', default=infinite
        '-e', '--evalue', 'set max evalue', default=0.01
        '-b', '--bscore', 'set min bit score', default=0
        '-qcov', '--qcov', 'set min qcov', default=0
        '-qcovh', '--qcovhsp', 'set min qcovhsp', default=0
        '-sort1', '--sort1', 'choose primary paramter to sort rows best to worst by entering string value: <qname|sname|pcntd|alnlen|msmtch|gpopen|qstart|q_end|sstart|s_end|evalue|bscore|qcov|qcovhsp>'
        '-sort2', '--sort2', 'choose secondary paramter to sort rows best to worst by entering string value: <qname|sname|pcntd|alnlen|msmtch|gpopen|qstart|q_end|sstart|s_end|evalue|bscore(default)|qcov|qcovhsp>'
        '-topq', '--gettopq', 'set number of top hits per query to output. Top hits based on preferred sort parameter - must use sort option'
    Library dependencies:
        'argparse', 'csv', 'operator'
        
