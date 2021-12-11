#!/usr/bin/env python3

import argparse

def GeneralArguments(defaults={
    'snakemake_executable': None,
    'tmpDir': None,
    'verbose':True,
}):
    """
    Arguments related to general settings
    """
    parser = argparse.ArgumentParser(add_help=False)
    general_args = parser.add_argument_group('general settings')
    general_args.add_argument("--snakemake_executable",
                         default=defaults.get('snakemake_executable'),
                         help="path to snakemake executable")
    general_args.add_argument("--tmpDir",
                         default=defaults.get('tmpDir'),
                         help="path to temporary directory")
    general_args.add_argument("-v", "--verbose",
                         dest="verbose",
                         action="store_true",
                         help="verbose output (default: '%(default)s')",
                         default=defaults["verbose"])

    return parser

