#!/usr/bin/env python3

__description__ = """
simple workflow for DNA mapping and coverage stats
usage example:
    ./pipeline.py -i /abs/path/to/input-dir/ -o output-dir hg19
"""

import argparse
import os
import sys
import textwrap
import parserCommon
import common_functions as cf

def parse_args(defaults={"verbose": False, "configFile": None,
                         "maxJobs": 5, "tmpDir": None, "trim": False,
                         "trimmer": "cutadapt", "trimmerOptions": "",
                         "dedup": True, "ext": ".fastq.gz",
                         "properPairs": True, "insertSizeMax": 1000,
                         "reads": ["_R1", "_R2"], "mapq": 0, 
                         "alignerOpts": "", "mateOrientation": "--fr",
                         "aligner":"Bowtie2"}):
    """
    Parse arguments from the command line.
    """
    GeneralArgs = parserCommon.GeneralArguments(defaults)
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__description__),
        parents=[GeneralArgs],
        add_help=True
    )

    # add positional argument specifying organism
    genomes = parserCommon.ListGenomes()
    parser.add_argument("genome", metavar="GENOME", 
    help="Genome acronym of the target organism. \
        Either a yaml file or one of: {}".format(", ".join(genomes)))
    # Workflow options

    # define required arguments 
    required = parser.add_argument_group('Required')

    required.add_argument("-i","--input-dir",
                          dest="indir",
                          help="specify path to fastq files",
                          required=True)

    required.add_argument("-o","--output-dir",
                          dest="outdir",
                          help="specify an output directory",
                          required=True)

    optional = parser.add_argument_group('Options')

    optional.add_argument("-c", "--configFile",
                         help="configuration file: config.yaml (default: '%(default)s')",
                         default=defaults["configFile"])

    optional.add_argument("-j", "--jobs",
                         dest="maxJobs",
                         metavar="INT",
                         help="maximum number of concurrently used cores (default: '%(default)s')",
                         type=int, default=defaults["maxJobs"])

    optional.add_argument("--alignerOpts",
                          help="Options that will be passed to Bowtie2 or bwa. You can specify things such as `--local` or "
                          "`--very-sensitive` here. The mate orientation and maximum insert size are specified "
                          "elsewhere. Read group information is set automatically. Note that you may need to escape "
                          "the first - (e.g., '\--very-fast'). Default: '%(default)s'.",
                          default=defaults["alignerOpts"])

    optional.add_argument("--mateOrientation",
                          help="The --fr, --ff, or --rf option for bowtie2 (default: '%(default)s')",
                          default=defaults["mateOrientation"])

    optional.add_argument("--dedup",
                          dest="dedup",
                          action="store_true",
                          help="retain only de-duplicated reads/read pairs "
                          "(given single-/paired-end data), recommended for "
                          "ChIP-seq data (default: '%(default)s')",
                          default=defaults["dedup"])

    optional.add_argument("--properPairs",
                          action="store_true",
                          help="retain only reads mapping in proper pairs (default: '%(default)s')",
                          default=defaults["properPairs"])

    optional.add_argument("--mapq",
                          dest="mapq",
                          metavar="INT",
                          help="retain only reads with at least the given "
                          "mapping quality. We recommend using"
                          "mapq of 3 or more for ChIP-seq to remove all true "
                          "multimapping reads. (default: '%(default)s')",
                          type=int,
                          default=defaults["mapq"])

    optional.add_argument("--insertSizeMax",
                          help="Maximum insert size allowed during mapping (default: '%(default)s')",
                          type=int,
                          default=defaults["insertSizeMax"])

    optional.add_argument("--aligner",
                          help="Program used for mapping: Bowtie2 or bwa (default: '%(default)s').",
                          choices=["Bowtie2","bwa"],
                          default=defaults["aligner"])

    return parser

def main():
    # set default values
    workflowDir, defaults = cf.setDefaults()

    # get command line arguments
    parser = parse_args(defaults)
    args = parser.parse_args()
    args, defaults = cf.handleUserArgs(args, defaults, parse_args)

    # create snakemake cmd that will be run as subprocess
    snakemake_cmd = cf.commonYAMLandLogs(workflowDir, defaults, args, __file__)

    # setup logging
    logfile_name = cf.logAndExport(args, os.path.basename(__file__))

    # run pipeline
    cf.runAndCleanup(args, snakemake_cmd, logfile_name)

if __name__ == "__main__":
    main()
    