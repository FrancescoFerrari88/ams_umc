# Snakemake pipeline for DNA mapping and coverage statistics

This pipeline allows to map a set of fastq files against a customizable reference genome and to compute the following coverage statistics:

- Mean per-base coverage over the whole reference genome;
- Median per-base coverage over the whole reference genome;
- Percentage of targeted base positions (PCT) in which coverage is greater than 30X;
- Mean coverage computed over exons, extended up- and down-stream by 6nt (extended exons);
- Aggregated gene-level mean coverage based on extended exons coverage;

## Installation

---

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
