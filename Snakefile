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
import subprocess
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
INDEX = config['INDEX']

# Full path to an uncompressed GFF file with known gene annotations.
GTF = config['GTF']
BED = config['BED']
CDS = config['CDS']

# Full path to a folder where final output files will be deposited.
OUT_DIR = config['OUT_DIR']
WORK_DIR = "/workdir"
HOME_DIR = config['HOME_DIR']

## set the working directory for each job (specific for CBSU qsub jobs)
USER = os.environ.get('USER')
JOB_ID = os.environ.get('JOB_ID')

# Samples and their corresponding filenames.
# paired-end:
FILES = json.load(open(config['SAMPLES_JSON'])) 
SAMPLES = sorted(FILES.keys())           
              
##--------------------------------------------------------------------------------------##
## RULES
##--------------------------------------------------------------------------------------##

## Create the final output directory if it doesn't already exist
if not os.path.exists(OUT_DIR):
            os.makedirs(OUT_DIR)

## Final expected output(s)
rule all: 
    input: 
      expand(join(OUT_DIR, 'CDS', '{sample}' + '.phase0.CDS.fasta'), sample = SAMPLES),
      expand(join(OUT_DIR, 'CDS', '{sample}' + '.phase1.CDS.fasta'), sample = SAMPLES),
      join(OUT_DIR, 'MultiQC', 'multiqc_report.html')
         


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
        # idx = rules.index.output.index,
        # bt2i = join(HOME_DIR, rstrip(DNA, '.fa') + '.ok')
    output: 
        bam = join(OUT_DIR, 'Bowtie2', '{sample}', '{sample}' + '.csorted.bowtie2.bam')
    log:
        join(OUT_DIR, 'Bowtie2', '{sample}', 'bowtie2.log')
    benchmark:
        join(OUT_DIR, 'Bowtie2', '{sample}', 'benchmark.tsv')
    message: 
        """--- Mapping PE sample "{wildcards.sample}" with Bowtie2."""
    run: 
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp {input.r1} {input.r2} ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp ' + join(INDEX, '*') + ' ' +  join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(DNA + '*') + ' ' +  join(WORK_DIR, USER, JOB_ID) +
                ' && cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && (bowtie2'                                     
                ' -p 16'   
                ' -x ' + os.path.basename(rstrip(DNA, '.fa')) +                    
                ' -1 {wildcards.sample}.R1.fq.gz' 
                ' -2 {wildcards.sample}.R2.fq.gz)> {log} 2>&1'
                ' | samtools sort -@ 8 -o csorted.bowtie2.bam -')
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, 'csorted.bowtie2.bam') + ' ' + join(OUT_DIR, 'Bowtie2', '{wildcards.sample}', '{wildcards.sample}' + '.csorted.bowtie2.bam'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

# Rule for phasing heterozygotes
rule phase_hets:
    input:
        dna = DNA,
        bam = join(OUT_DIR, 'Bowtie2', '{sample}', '{sample}' + '.csorted.bowtie2.bam')
    output:
        phase1 = join(OUT_DIR, 'Bowtie2', '{sample}', '{sample}' + '.phased.0.bam'),
        phase2 = join(OUT_DIR, 'Bowtie2', '{sample}', '{sample}' + '.phased.1.bam'),
        chimer = join(OUT_DIR, 'Bowtie2', '{sample}', '{sample}' + '.phased.chimera.bam')
    log:
        join(OUT_DIR, 'Bowtie2', 'logs', '{sample}' + '.phase.log')
    benchmark:
        join(OUT_DIR, 'Bowtie2', 'logs', '{sample}' + '.phase.benchmark.tsv')
    message: 
        """--- Phasing heterozygotes for sample "{wildcards.sample}". """
    run:
        # Extract a sequence for each transcript in the GTF file.
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp ' + join(INDEX, '*') + ' ' +  join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(DNA + '*') + ' ' +  join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.bam} ' + join(WORK_DIR, USER, JOB_ID))
        shell('cd ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && samtools phase'
                ' -b {wildcards.sample}.phased'
                ' -Q 20'
                ' -D 200'
                ' --reference ' + os.path.basename(DNA) +
                ' {wildcards.sample}.csorted.bowtie2.bam > {log} 2>&1')  
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.phased*') + ' ' + join(OUT_DIR, 'Bowtie2', '{wildcards.sample}'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))


# Rule for generating pileup and extracting consensus fasta
rule mpileup_fasta_phase0:
    input:
        dna = DNA,
        bam = join(OUT_DIR, 'Bowtie2', '{sample}', '{sample}' + '.phased.0.bam')
    output:
        fasta = join(OUT_DIR, 'CDS', '{sample}' + '.phase0.CDS.fasta')
    params:
        cds = CDS,
        bed = BED,
        gtf = GTF
    log:
        join(OUT_DIR, 'CDS', 'logs', '{sample}' + '.phase0.log')
    benchmark:
        join(OUT_DIR, 'CDS', 'logs', '{sample}' + '.phase0.benchmark.tsv')
    message: 
        "--- Generating pileup and extracting phase0 genome consensus fasta for sample {wildcards.sample}"
    run:
        # Extract a sequence for each transcript in the GTF file.
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp ' + join(INDEX, '*') + ' ' +  join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(DNA + '*') + ' ' +  join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.bam} {params.cds} {params.bed} {params.gtf} ' + join(WORK_DIR, USER, JOB_ID))
        transcripts = [line.strip('\n') for line in open(CDS)]
        for line in transcripts:
            try:
                shell('cd ' + join(WORK_DIR, USER, JOB_ID) + ' && grep {line} ' + os.path.basename(BED) + ' | grep CDS | samtools mpileup -u -l - -f ' + os.path.basename(DNA) + ' {wildcards.sample}.phased.0.bam'
                    ' | bcftools call -c | vcfutils.pl vcf2fq | seqtk seq -A /dev/stdin/ > {line}.fa' 
                    ' && grep {line} ' + os.path.basename(GTF) + ' | grep CDS | gffread - -g {line}.fa -x {line}.CDS.fa'
                    ' && cat {line}.CDS.fa >> ' +  join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.phase0.CDS.fasta') +
                    ' && rm {line}.fa.fai {line}.fa {line}.CDS.fa')
            except Exception:
                continue
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.phase0.CDS.fasta') + ' ' + join(OUT_DIR, 'CDS'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))

# Rule for generating pileup and extracting consensus fasta
rule mpileup_fasta_phase1:
    input:
        dna = DNA,
        bam = join(OUT_DIR, 'Bowtie2', '{sample}', '{sample}' + '.phased.1.bam')
    output:
        fasta = join(OUT_DIR, 'CDS', '{sample}' + '.phase1.CDS.fasta')
    params:
        cds = CDS,
        bed = BED,
        gtf = GTF
    log:
        join(OUT_DIR, 'CDS', 'logs', '{sample}' + '.phase1.log')
    benchmark:
        join(OUT_DIR, 'CDS', 'logs', '{sample}' + '.phase1.benchmark.tsv')
    message: 
        "--- Generating pileup and extracting phase1 genome consensus fasta for sample {wildcards.sample}"
    run:
        # Extract a sequence for each transcript in the GTF file.
        shell('mkdir -p ' + join(WORK_DIR, USER, JOB_ID) + 
                ' && cp ' + join(INDEX, '*') + ' ' +  join(WORK_DIR, USER, JOB_ID) +
                ' && cp ' + join(DNA + '*') + ' ' +  join(WORK_DIR, USER, JOB_ID) +
                ' && cp {input.bam} {params.cds} {params.bed} {params.gtf} ' + join(WORK_DIR, USER, JOB_ID))
        transcripts = [line.strip('\n') for line in open(CDS)]
        for line in transcripts:
            try:
                shell('cd ' + join(WORK_DIR, USER, JOB_ID) + ' && grep {line} ' + os.path.basename(BED) + ' | grep CDS | samtools mpileup -u -l - -f ' + os.path.basename(DNA) + ' {wildcards.sample}.phased.1.bam'
                    ' | bcftools call -c | vcfutils.pl vcf2fq | seqtk seq -A /dev/stdin/ > {line}.fa' 
                    ' && grep {line} ' + os.path.basename(GTF) + ' | grep CDS | gffread - -g {line}.fa -x {line}.CDS.fa'
                    ' && cat {line}.CDS.fa >> ' +  join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.phase1.CDS.fasta') +
                    ' && rm {line}.fa.fai {line}.fa {line}.CDS.fa')
            except Exception:
                continue
        shell('mv ' + join(WORK_DIR, USER, JOB_ID, '{wildcards.sample}.phase1.CDS.fasta') + ' ' + join(OUT_DIR, 'CDS'))
        shell('rm -r ' + join(WORK_DIR, USER, JOB_ID))


## Rule to collate fastQC and Bowtie2 outputs with multiQC
rule multiQC:
    input:
        expand(join(OUT_DIR, 'Bowtie2', '{sample}', '{sample}' + '.csorted.bowtie2.bam'), sample = SAMPLES),
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