# Snakemake pipeline for DNA mapping and coverage statistics

This pipeline allows to map a set of fastq files against a customizable reference genome and to compute the following coverage statistics:

- Mean per-base coverage over the whole reference genome;
- Median per-base coverage over the whole reference genome;
- Percentage of targeted base positions (PCT) in which coverage is greater than 30X;
- Mean coverage computed over exons, extended up- and down-stream by 6nt (extended exons);
- Aggregated gene-level mean coverage based on extended exons coverage;

## Installation

To use this tool, you need git and conda.

1.  First, clone this repository:

        git clone git@github.com:FrancescoFerrari88/ams_umc.git

2.  from inside the cloned repository, create a conda environment using the environment.yaml file

        conda env create -f environment.yaml

3.  Next, you can activate your environment using:

        conda activate snakemake_amsumc

You can test whether the creation of the conda environment was successfull by checking the output of

        conda env list

and

        ./pipeline.py -h

## Generate Supporting Organism Files

Before using the pipeline, you need to create an organim folder with the genome of interest, Bowtie2 index of the same genome, a gtf gene annotation and a bed file with non-overlapping extended exons (6nt up- and down-stream).

### Example

0. Create a folder called GRCh37p13
   mkdir /path/to/GRCh37p13; cd /path/to/GRCh37p13; mkdir genome bowtie_index annotation
1. In /path/to/GRCh37p13/genome, download NCBI human genome GRCh37.p13 (hg19) that you find at the following address:
   > https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/

After downloading all chromosomes, concatenate them to generate a genome fasta file

        cat \*.gz > GRCh37.p13.fa.gz

In the same folder, generate an unzipped version of the genome with its index:

        mkdir fasta; cd fasta; gunzip ../GRCh37.p13.fa.gz; samtools faidx GRCh37.p13.fa

2.  Generate Bowtie2-index for the genome GRCh37.p13. To this end, go into the folder bowtie_index that you previously created, and run

        bowtie2-build --threads 6 ../genome/GRCh37.p13.fa.gz GRCh37.p13

This will generate the index. The process may take some time.

3.  Download a gtf gene model annotation compatible with your genome (name of chromosomes must match!). To download the corresponding gtf annotation for NCBI GRCh37.p13, go into the folder annotation, and run

        wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gtf.gz

Next, sort and index the gtf file (e.g using IGV tools).
Finally, you can create a bed file with non-overlapping exons extended by 6 nucleotides up- and down-stream for the downloaded annotation. For the annotation at hand, you can do that by runnign the following command:

        awk '{if ($3=="exon"){print $0}}' GRCh37*latest_genomic.sorted.gtf | awk 'BEGIN {OFS="\t"} {print $1,$4-6,$5+6,$10,$6,$7 }' | grep "^NC*" | bedtools merge -c 4 -o distinct | awk '! /;,/{print $0}' | awk -v OFS="\t" '{print $1,$2,$3,substr($4,2,length($4)-3) }' > uniq_merged_exons.bed

## Organim Config

Next, you need to configure the organims/<organim>.yaml file to point to the respective files/folders.

1. Customize organism config file:

inside the repository, open organisms/<organim>.yaml and fill in the absolute paths to fasta_genome (1), bowtie2_index (2), extended_exon_bed (3), genome_length (4). To get the total genome length, run:

        bowtie2-inspect -s /path/to/bwt2-idx/bt2_base | awk '{ FS = "\t" } ; BEGIN{L=0}; {L=L+$3}; END{print L}'

## Default Config
