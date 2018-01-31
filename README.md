tools and pipelines for deep sequencing analysis

# Mutant Discovery
This pipeline was optimised for mutagen induced SNP discovery in sequence data from exome capture experiments. However, flexible options allow wide application for variant analysis in any deep sequencing data.

## Tool descriptions

**Noisefinder.pyc**
Regions rich in mismatches/poor coverage after read alignment can often signify misalignment or mixed alignment due to allelism, polyploidy, or presence/absemce variation. Given a pileup file, noisefinder reports regions containing a density of SNVs above a user defined threshold over a given min length and min read depth (prints to STDOUT).
    
**SNPlogger.pyc**
Will parse an mpileup file and log all SNPs and indels that satisfy parameters. Final tally printed to STDOUT. Outfile is formatted as tab sep fields: <seqid> <position(1based)> <polymorphic-type> <frequency(float)>. Compatible with STDIN.

**SNPtracker.pyc**
Finds sequence IDs/regions with coinciding polymorphic features across multiple SNPlogger generated files

## Example Workflow
This specific workflow is designed to discover sequences/contigs that contain mutagen induced variation occuring independently across a number of mutants. In mutagenesis experiments for which single gene knockouts can be selected for phenotypically, such a finding is strongly indidcative that the target gene has been isolated given a sufficient number of mutants. It is based on generating a de novo assembly from wild-type NGS reads and then aligning mutant NGS reads independently against the wild-type assembly and recording any mismatches between each mutant and the wild-type. Ideally, the wild-type should be parental to the mutants and all be near-isogenic lines in order to minimise noise due to normal genetic variation. This pipeline is inspired by similar pipelines such as MutantHunter (https://github.com/steuernb/MutantHunter), but with an alternate approach and added flexibility.

### Prerequisites
**Python 2.7**
see https://www.python.org/download/releases/2.7/
**BWA or other suitable aligner**
see http://bio-bwa.sourceforge.net/
**Samtools 1.5 or later**
see http://samtools.sourceforge.net/
**_De novo_ assembly software**
can be of your choosing but should be relatively stringent to avoid collapsing highly homologous yet distinct sequences as a single consensus. The CLC assembly cell (https://www.qiagenbioinformatics.com/products/clc-assembly-cell/) or MaSuRCA assembler (http://www.genome.umd.edu/masurca.html) generate suitable assemblies.

### Preprocessing

#### clean raw data
if data is not already cleaned/trimmed, software such as Trimmomatic can be used (http://www.usadellab.org/cms/?page=trimmomatic)
for example if your data is paired-end then for each fastq file you would do something like:

	trimmomatic PE -threads 8 -phred33 readsIn_1.fq.gz readsIn_2.fq.gz readsOut_1.clean.fq.gz readsOut_1.unpaired.fq.gz readsOut_2.clean.fq.gz readsOut_2.unpaired.fq.gz ILLUMINACLIP:adapter_seqs.fasta:2:30:10:8:TRUE LEADING:28 TRAILING:28 MINLEN:20

#### _de novo_ assembly of wild-type 
as stated earlier, use tool of your choice, but ensure a decent N50 as a quality assembly is the crux of the whole procedure. Note that some assemblers prefer the input to be raw, uncleaned data (MaSuRCA).

#### mapping to wild-type
as well as te mutants, the wild-type reads should also be mapped to their assembly to ensure greater accuracy in later steps. Again, choice of alignment software is at your discretion but BWA and samtools offer a straightforward process. Initially index the assembly:

	bwa index WT_assembly.fasta
	samtools faidx WT_assembly.fasta

then run the following steps for each mutant and wild-type:

	bwa aln assembly.fasta read1.fastq > read1.aln
	bwa aln assembly.fasta read2.fastq > read2.aln
	bwa sampe assembly.fasta read1.aln read2.aln read1.fastq read2.fastq > raw.sam
	samtools view -f2 -Shub -o raw.bam raw.sam
	samtools sort raw.bam sorted
	samtools rmdup sorted.bam rmdup.bam
	samtools index rmdup.bam

#### Mutant discovery steps
Using example files for wildtype "WT.rmdup.bam" plus mutants "mut1.rmdup.bam", "mut2.rmdup.bam", "mut3.rmdup.bam"...

**1) create pileup of WT bam only.**

    samtools mpileup -BQ0 -aa -f WT_assembly.fasta WT.rmdup.bam > WT.pileup
       
**2) run Noisefinder on WT pileup.**

    python Noisefinder.pyc -i WT.pileup > WT.noise.log

**3) run SNPlogger on WT and mutants using WT.noise.log to mask rubbish regions.**
Run "python SNPlogger.pyc -h" to see additional options as noise.log files generated from mutants can be used as features themselves in later steps.

for WT, SNPlogger can be run on already created pileup:

    python SNPlogger.pyc -i WT.pileup -b WT.noise.log -o WT.snp.log

for mutants, pileups can be created on the fly and piped directly to SNPlogger.

one by one:

    samtools mpileup -aa -BQ0 -f WT_assembly.fasta mut1.rmdup.bam | python SNPlogger.pyc -b WT.noise.log -o mut1.snp.log
or in a loop:

    for i in mut{1..3}; do samtools mpileup -aa -BQ0 -f WT_assembly.fasta ${i}.rmdup.bam | python SNPlogger.pyc -b WT.noise.log -o ${i}.snp.log; done

note that SNPlogger will print a summary of SNP statistics to the screen upon completion of each file. To save these stats, use ">" to redirect standard output to a file:

    for i in m{1..3}; do samtools mpileup -aa -BQ0 -f WT_assembly.fasta ${i}.bam | python SNPlogger.pyc -b WT.noise.log -o ${i}.snp.log > ${i}.stats.txt; done

**4) run SNPtracker on snp.log files.**
used "-w" for WT file(s), "-m" for mutant files. SNPtracker can still work without a WT or with >1 WT. This step is relatively fast and can complete in seconds:

    python SNPtracker.pyc -w WT.snp.log -m mut1.snp.log mut2.snp.log mut3.snp.log

run "python SNPtracker.pyc -h" to see options. E.g. to only report polymorphisms that are C>T, G>A or indels, use "-s" option. To filter polymorphisms based on frequency (default=0.8), use "-f" option. To create detailed reports in addition to the summary report, use "-v T" option. To tolerate N mutants with identical SNPs (e.g. siblings), use "-t N" option:

    python SNPtracker.pyc -w WT.snp.log -m mut1.snp.log mut2.snp.log mut3.snp.log -s C\>T G\>A indel -f 0.9 -v T -t 2

This will create a "SNPtracker.summary" report whose output looks like the following:
	
>	# wildtypes: WT.snp.log
>	# mutants: mut1.snp.log, mut2.snp.log, mut3.snp.log
>	# selected: C>T, G>A, indel
>	# filtered: 0.9
>	# proximal: OFF
>	# tolerate: 2
>	
>	### polymorphic in 3 mutants ###
>	contig_5565     (mut1.snp.log, mut2.snp.log, mut3.snp.log)
>	
>	### polymorphic in 2 mutants ###
>	contig_8725     (mut1.snp.log, mut3.snp.log)
>	contig_1252     (mut2.snp.log, mut3.snp.log)

If the "-v" option is set to true, the detailed reports generated display the coordinate and type of polymorphisms found in each of the mutants for a given candidate contig/sequence. A report is generated for each N>1. i.e. if there are 5 mutants, SNPtracker will create an N5.report, N4.report, N3.report and N2.report. That is unless a lower limit is specified with the "-n" option. An example of an N2.report may look like the following:

>	### poymorphic in 2 mutants ###
>
>	<contig_8725>
>	mut1.snp.log [(1200, G>A, 0.98)]
>	mut3.snp.log [(3210, C>T, 0.91), (4128, indel+5, 1.0)]
>
>	<contig_1252>
>	mut2.snp.log [(852, C>T, 0.95), (3069, G>A. 0.90)]
>	mut3.snp.log [(2567, indel-8, 1.0)]

Alignments for these candidate contigs can also be inspected visually upon loading the bam files into a genome browser such as IGV (http://software.broadinstitute.org/software/igv/).

SNPtracker also provides a "-p/--proximal" option that tells SNPtracker to find features that reside close to each other within a user defined window (min=1000 bases) rather than only coinciding on a particular contig. This is suitable if working with large scaffolds or assemblies from long read sequencing such as PacBio.

# Miscellaneous

        
