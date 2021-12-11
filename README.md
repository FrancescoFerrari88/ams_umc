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

0.  Create a folder called GRCh37p13

        mkdir /path/to/GRCh37p13; cd /path/to/GRCh37p13; mkdir genome bowtie_index annotation

1.  In /path/to/GRCh37p13/genome, download NCBI human genome GRCh37.p13 (hg19) that you find at the following address:

    > https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/

After downloading all chromosomes, concatenate them to generate a genome fasta file

        cat *.gz > GRCh37.p13.fa.gz

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

## Organism Config

Next, you need to configure the organisms/organism.yaml file to point to the respective files/folders.

1. Customize organism config file:

Inside the repository, open organisms/organism.yaml and fill in the absolute paths to fasta_genome (1 - non gzipped), bowtie2_index (2), extended_exon_bed (3), genome_length (4). To get the total genome length, run:

        bowtie2-inspect -s /path/to/bwt2-idx/bt2_base | awk '{ FS = "\t" } ; BEGIN{L=0}; {L=L+$3}; END{print L}'

An example organim config file is present in organisms/hg19.yaml

## Usage

To have an overview of the usage and available parameters, run

        ./pipeline.py -h

        usage: ./pipeline.py [-h] [--snakemake_executable SNAKEMAKE_EXECUTABLE] [--tmpDir TMPDIR] [-v] -i INDIR -o OUTDIR
                     [-c CONFIGFILE] [-j INT] [--alignerOpts ALIGNEROPTS] [--mateOrientation MATEORIENTATION] [--dedup]
                     [--properPairs] [--mapq INT] [--insertSizeMax INSERTSIZEMAX] [--aligner {Bowtie2,bwa}]
                     GENOME

        simple workflow for DNA mapping and coverage stats
        usage example:
            ./pipeline.py -i /abs/path/to/input-dir/ -o output-dir hg19

        positional arguments:
        GENOME                Genome acronym of the target organism. Either a yaml file or one of: hg19

        optional arguments:
        -h, --help            show this help message and exit

        general settings:
        --snakemake_executable SNAKEMAKE_EXECUTABLE
                                path to snakemake executable
        --tmpDir TMPDIR       path to temporary directory
        -v, --verbose         verbose output (default: 'True')

        Required:
        -i INDIR, --input-dir INDIR
                                specify path to fastq files
        -o OUTDIR, --output-dir OUTDIR
                                specify an output directory

        Options:
        -c CONFIGFILE, --configFile CONFIGFILE
                                configuration file: config.yaml (default: 'None')
        -j INT, --jobs INT    maximum number of concurrently used cores (default: '5')
        --alignerOpts ALIGNEROPTS
                                Options that will be passed to Bowtie2 or bwa. You can specify things such as `--local` or `--very-
                                sensitive` here. The mate orientation and maximum insert size are specified elsewhere. Read group
                                information is set automatically. Note that you may need to escape the first - (e.g., '\--very-fast').
                                Default: 'None'.
        --mateOrientation MATEORIENTATION
                                The --fr, --ff, or --rf option for bowtie2 (default: '--fr')
        --dedup               retain only de-duplicated reads/read pairs (given single-/paired-end data), recommended for ChIP-seq
                                data (default: 'False')
        --properPairs         retain only reads mapping in proper pairs (default: 'False')
        --mapq INT            retain only reads with at least the given mapping quality. We recommend usingmapq of 3 or more for ChIP-
                                seq to remove all true multimapping reads. (default: '0')
        --insertSizeMax INSERTSIZEMAX
                                Maximum insert size allowed during mapping (default: '1000')
        --aligner {Bowtie2,bwa}
                                Program used for mapping: Bowtie2 or bwa (default: 'Bowtie2').

## Pipeline Parameters

You can control the behaviour of the pipeline, by changing the values of the parameters. You can do that in 3 ways:

1. Globally, by changing the values of the default config file in config/defaults.yaml;
2. Per single run, providing a custom config file using the option --configFile;
3. Per single option, providing value to single option flags;

Some parameters can only be controlled using the config file (either default or custom). Of particular importance is the parameter "reads".

> "reads" specify the format of the portion of the fastq filename that carries information as to whether the read is R1 or R2 in a paired-end experiment. This portion of the filename is expected to be immediately before the extension (ext: "fastq.gz") and immediately after the sample name.
>
> ### Example
>
> If the paired fastq file names are "D1_S1_L001_R1_001.fastq.gz" (R1) and "D1_S1_L001_R2_001.fastq.gz" (R2), the parameter "reads" will be:
>
>        reads: ['R1_001','R2_001']

## Understanding the Output

The pipeline generates 5 folders in the output path:

- **originalFASTQ**: collects soft links to the original input fastq files
- **Bowtie2**: collects summary statistics about the mapping of the input samples
- **CRAM**: collects alignment data in cram format
- **CoverageStats**: collects per-base coverage stats computed over the whole genome
- **ExonGeneCoverage**: collects exon mean coverage ( _sample_name_\_exonCov.txt) and gene-level aggregate mean coverage (_sample_name_\_aggGeneCov.txt)

In the output folder, you will also find copy of the config parameters used in the analysis (pipeline.config.yaml) and the complete stdout generated during pipeline execution (pipeline.py_run-1.yaml). If you run the pipeline several times, new stdout log files will be added.
