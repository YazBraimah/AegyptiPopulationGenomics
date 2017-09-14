"""
Author: Y. Ahmed-Braimah
--- Snakemake workflow for the assembly and CDS extraction 
--- of Ae. aegypti exome data
"""

import json
import os
from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output

##--------------------------------------------------------------------------------------##
## Functions
##--------------------------------------------------------------------------------------##

# To print process messages
def message(x):
  print()

# To remove suffix from a string
def rstrip(text, suffix):
    if not text.endswith(suffix):
        return text
    return text[:len(text)-len(suffix)]

## define environment variables

##--------------------------------------------------------------------------------------##
## Global config files: 
##--------------------------------------------------------------------------------------##

configfile: 'config.yml'

# Full path to an uncompressed FASTA file with all chromosome sequences.
DNA = config['DNA']

# Full path to an uncompressed GFF file with known gene annotations.
GTF = config['GTF']

# Full path to a folder where final output files will be deposited.
OUT_DIR = config['OUT_DIR']

## set the working directory for each job (specific for CBSU qsub jobs)
USER = os.environ.get('USER')
JOB_ID = os.environ.get('JOB_ID')
WORK_DIR = "/workdir"

HOME_DIR = config['HOME_DIR']
# message("The current working directory is " + WORK_DIR)

# Samples and their corresponding filenames.
# paired-end:
FILES = json.load(open(config['SAMPLES_JSON'])) 
SAMPLES = sorted(FILES.keys())           
              
##--------------------------------------------------------------------------------------##
## RULES
##--------------------------------------------------------------------------------------##

## Final expected output(s)
rule all: 
    input: 
      expand(join(OUT_DIR, 'Bowtie2', '{sample}',  '.csorted.bowtie2.bam'), sample = SAMPLES),
      expand(join(OUT_DIR, 'fastQC', '{sample}', '{sample}' + '.R2_fastqc.html'), sample = SAMPLES),
      join(OUT_DIR, 'MultiQC', 'multiqc_report.html')
         

## Rule to generate bowtie2 genome index 
rule index:
    input:
        dna = DNA
    output:
        index = join(dirname(DNA), rstrip(DNA, '.fa') + '.rev.1.bt2'),
        bt2i = join(dirname(DNA), rstrip(DNA, '.fa') + '.ok')
    log:
        join(dirname(DNA), 'bt2.index.log')
    benchmark:
        join(dirname(DNA), 'bt2.index.benchmark.tsv')
    message: 
        """--- Building bowtie2 genome index """
    run:
        shell('samtools faidx {input.dna}')
        shell('bowtie2-build {input.dna} ' + join(dirname(DNA), rstrip(DNA, '.fa')) + ' > {log} 2>&1')
        shell('touch ' + join(dirname(DNA), rstrip(DNA, '.fa') + '.ok'))


# Rule to check PE read quality
rule fastqc:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output:
        r1 = join(OUT_DIR, 'fastQC', '{sample}', '{sample}' + '.R1_fastqc.html'),
        r2 = join(OUT_DIR, 'fastQC', '{sample}', '{sample}' + '.R2_fastqc.html')
    log:
        join(OUT_DIR, 'fastQC', '{sample}', 'fastQC_pe.log')
    benchmark:
        join(OUT_DIR, 'fastQC', '{sample}', 'fastQC_pe.benchmark.tsv')
    message: 
        """--- Checking trimmed read quality of sample "{wildcards.sample}" with FastQC """
    run:
        if not os.path.exists(join(OUT_DIR, 'fastQC')):
            os.makedirs(join(OUT_DIR, 'fastQC'))

        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) +
              ' && cp {input.r1} {input.r2} ' + join(WORK_DIR, USER, JOB_ID) +
              ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && fastqc {wildcards.sample}.R1.fq.gz {wildcards.sample}.R2.fq.gz' 
                ' > {log} 2>&1')
        shell('cd ' + join(WORK_DIR, USER, JOB_ID) + ' && rm {wildcards.sample}.R1.fq.gz {wildcards.sample}.R2.fq.gz')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID) + '/* ' + join(OUT_DIR, 'fastQC', '{wildcards.sample}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))
        

## Rule for mapping PE reads to the genome with Bowtie2
rule Bowtie2:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2'],
        idx = rules.index.output.bt2i
    output: 
        bam = join(OUT_DIR, 'Bowtie2', '{sample}', 'csorted.bowtie2.bam')
    params: 
        gtf = GTF
    log:
        join(OUT_DIR, 'Bowtie2', '{sample}', 'bowtie2.log')
    benchmark:
        join(OUT_DIR, 'Bowtie2', '{sample}', 'benchmark.tsv')
    message: 
        """--- Mapping PE sample "{wildcards.sample}" with Bowtie2."""
    run: 
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp {input.r1} {input.r2} ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp ' + join(dirname(DNA), rstrip(DNA, '.fa') + '*') + ' ' +  join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && bowtie2'                                     
                ' -p 16'   
                ' -x ' + os.path.basename(join(dirname(DNA), rstrip(DNA, '.fa'))) +                    
                ' -1 {wildcards.sample}.R1.fq.gz' 
                ' -2 {wildcards.sample}.R2.fq.gz'
                ' | samtools view -bS - > bowtie2.bam'
                ' > {log} 2>&1')
        shell('samtools sort bowtie2.bam -o csorted.bowtie2.bam')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, 'csorted.bowtie2.bam') + ' ' + join(OUT_DIR, 'Bowtie2', '{sample}', 'csorted.bowtie2.bam')
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))


# Rule for generating pileup and extracting consensus fasta
rule mpileup_fasta:
    input:
        dna = DNA,
        bam = rules.Bowtie2.output.bam
    output:
        fasta = join(OUT_DIR, 'genome_fasta', '{sample}' + '.fasta')
    log:
        pileup = join(OUT_DIR, 'genome_fasta', 'logs', '{sample}' + '.log')
    benchmark:
        join(OUT_DIR, 'genome_fasta', 'logs', '{sample}' + '.log')
    message: 
        "--- Generating pileup and extracting genome consensus fasta"
    run:
        # Extract a sequence for each transcript in the GTF file.
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp {input.bam} ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp ' + join(dirname(DNA), rstrip(DNA, '.fa') + '*') + ' ' +  join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && samtools mpileup'
                ' -uf ' + os.path.basename(join(dirname(DNA), rstrip(DNA, '.fa'))) + 
                ' {input.bam}'
                ' | bcftools call -c | vcfutils.pl vcf2fq > {wildcards.sample}.fq')
        shell('seqtk seq -A {wildcards.sample}.fq > {wildcards.sample}.fasta')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.fasta') + ' ' + join(OUT_DIR, 'genome_fasta')
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))


## Rule to collate fastQC and Bowtie2 outputs with multiQC
rule multiQC:
    input:
        expand(join(OUT_DIR, 'Bowtie2', '{sample}', 'csorted.bowtie2.bam'), sample = SAMPLES),
        expand(join(OUT_DIR, 'fastQC', '{sample}', '{sample}' + '.R2_fastqc.html'), sample = SAMPLES)
    output:
        file = join(OUT_DIR, 'MultiQC', 'multiqc_report.html')
    log:
        join(OUT_DIR, 'MultiQC', 'multiQC.log')
    benchmark:
        join(OUT_DIR, 'MultiQC', 'multiQC.benchmark.tsv')
    message: 
        """--- Running MultiQC """
    run:
        shell('ls -1 ' + join(OUT_DIR) + '/Bowtie2/*/bowtie2.log > ' + join(OUT_DIR, 'summary_files.txt'))
        shell('ls -1 ' + join(OUT_DIR) + '/fastQC/*/*fastqc.zip >> ' + join(OUT_DIR, 'summary_files.txt'))
        shell('multiqc'
                ' -f'
                ' -o ' + join(OUT_DIR, 'MultiQC') + ' -d -dd 2 -l ' + join(OUT_DIR, 'summary_files.txt') +
                ' > {log} 2>&1')

## Rule to extract CDS sequences from genome fasta
rule gffread:
	input:
		fasta = join(OUT_DIR, 'genome_fasta', '{sample}' + 'fasta')
	params:
		gtf = GTF
	log:
        join(OUT_DIR, 'Bowtie2', '{sample}', 'bowtie2.log')
    benchmark:
        join(OUT_DIR, 'Bowtie2', '{sample}', 'benchmark.tsv')
    message: 
        """--- Mapping PE sample "{wildcards.sample}" with Bowtie2."""
    run: 