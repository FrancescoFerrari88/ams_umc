#!/usr/bin/env python

# general functions for workflow #############################################
##############################################################################
import subprocess
import os
import yaml
import re
import sys
import glob

def merge_dicts(x, y):
    z = {}
    z = x.copy()
    if y:
        z.update(y)
    return z

def load_configfile(configFiles, verbose, info='Config'):
    with open(configFiles, "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    config = config

    if verbose:
        print("\n--- " + info + " ---------------------------------------------------------------------")
        print("config file: {}".format(configFiles))
        for k, v in sorted(config.items()):
            print("{}: {}".format(k, v))
        print("-" * 80, "\n")
    return config

def setDefaults():
    """
    Set a number of variables used in the wrappers and the defaults
    """
    # workflowDir is located where this current file is
    workflowDir = os.path.dirname(__file__)

    # defaults
    defaults = load_configfile(os.path.join(workflowDir, "config", "defaults.yaml"), False)

    return workflowDir, defaults

def handleUserArgs(args, defaults, args_func):
    """
    If a user supplies a custom YAML file then that must replace the defaults.
    However command line options need to take precedence, so update defaults and
    simply reparse the command line options (with args_func().parse_args())
    """
    if args.configFile:
        if not os.path.exists(args.configFile):
            sys.exit("\nError! Provided configFile (-c) not found! ({})\n".format(args.configFile))
        user_config = load_configfile(args.configFile, False)
        defaults = merge_dicts(defaults, user_config)
        parser = args_func(defaults)
        args = parser.parse_args()
    defaults.update(vars(args)) # here vars returns the __dict__ method of args
    return args, defaults

def write_configfile(configFile, config):
    with open(configFile, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

def commonYAMLandLogs(workflowDir, defaults, args, callingScript):
    """
    Merge dictionaries, write YAML files, construct the snakemake command
    and create the DAG
    """
    workflowName = ".".join(os.path.basename(callingScript).split(".")[:-1])

    # create outdit and temporary directory 
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.tmpDir, exist_ok=True)

    # save to configs.yaml in outdir
    config = defaults
    config.update(vars(args))  # This modify args after handling a user config file to still make it to the YAML given to snakemake!
    write_configfile(os.path.join(args.outdir, '{}.config.yaml'.format(workflowName)), config)

    snakemake_cmd = """
                    TMPDIR={tempDir} PYTHONNOUSERSITE=True {snakemake} {snakemakeOptions} \
                    --latency-wait {latency_wait} --snakefile {snakefile} --jobs {maxJobs} \
                    --directory {workingdir} --configfile {configFile} --keep-going
                    """.format(snakemake=args.snakemake_executable,
                               latency_wait=10,
                               snakefile=os.path.join(workflowDir, "Snakefile"),
                               maxJobs=args.maxJobs,
                               workingdir=args.outdir,
                               snakemakeOptions='',
                               tempDir=args.tmpDir,
                               configFile=os.path.join(args.outdir, '{}.config.yaml'.format(workflowName))).split()

    if args.verbose:
        snakemake_cmd.append("--printshellcmds")

    return " ".join(snakemake_cmd)

def logAndExport(args, workflowName):
    """
    Set up logging
    """
    # Write snakemake_cmd to log file
    fnames = glob.glob(os.path.join(args.outdir, '{}_run-[0-9]*.log'.format(workflowName)))
    if len(fnames) == 0:
        n = 1  # no matching files, this is the first run
    else:
        fnames.sort(key=os.path.getctime)
        n = int(fnames[-1].split("-")[-1].split(".")[0]) + 1  # get new run number
    # append the new run number to the file name
    logfile_name = "{}_run-{}.log".format(workflowName, n)

    return logfile_name


def runAndCleanup(args, cmd, logfile_name):
    """
    Actually run snakemake. Kill its child processes on error.
    Also clean up when finished.
    """
    if args.verbose:
        print("\n{}\n".format(cmd))

    # write log file
    f = open(os.path.join(args.outdir, logfile_name), "w")
    f.write(" ".join(sys.argv) + "\n\n")
    f.write(cmd + "\n\n")

    # Run snakemake, stderr -> stdout is needed so readline() doesn't block
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if args.verbose:
        print("PID:", p.pid, "\n")

    while p.poll() is None:
        stdout = p.stdout.readline(1024)
        if stdout:
            sys.stdout.write(stdout.decode('utf-8'))
            f.write(stdout.decode('utf-8'))
            sys.stdout.flush()
            f.flush()

    # This avoids the race condition of p.poll() exiting before we get all the output
    stdout = p.stdout.read()
    if stdout:
        sys.stdout.write(stdout.decode('utf-8'))
        f.write(stdout.decode('utf-8'))
    f.close()

    # Exit with an error if snakemake encountered an error
    if p.returncode != 0:
        sys.stderr.write("Error: snakemake returned an error code of {}, so processing is incomplete!\n".format(p.returncode))
    else:
        if os.path.exists(os.path.join(args.outdir, ".snakemake")):
            import shutil
            shutil.rmtree(os.path.join(args.outdir, ".snakemake"), ignore_errors=True)

def load_organism_data(genome, maindir, verbose):

    if os.path.isfile(os.path.join(maindir, 'organisms', genome + ".yaml")):
        organism = load_configfile(os.path.join(maindir, 'organisms', genome + ".yaml"), verbose, "Genome")
    else:
        exit("ERROR: Genome configuration file NOT found for: {}\n".format(genome))
    return organism

def get_sample_names(infiles, ext, reads):
    """
    Get sample names without file extensions
    """
    s = set()
    lext = len(ext)
    l0 = len(reads[0])
    l1 = len(reads[1])
    for x in infiles:
        x = os.path.basename(x)[:-lext]
        if x.endswith(reads[0]):
            x = x[:-l0]
            s.add(x)
        elif x.endswith(reads[1]):
            x = x[:-l1]
            s.add(x)
        else:
            sys.stderr.write("Warning! {} does not have {} as its name suffix. "
                             "Either change it or modify the 'reads' in the "
                             "config.yaml to your deired ones.\n".format(x, reads))

    if sorted(list(s)) == []:
        sys.exit("Error! No sample has the right read suffix ({}). "
                 "Please modify them or update the config.yaml with "
                 "your desired suffix.".format(reads))
    return sorted(list(s))

def is_paired(infiles, ext, reads):
    """
    Check for paired-end input files
    """
    pairedEnd = False
    infiles_dic = {}
    for infile in infiles:
        fname = os.path.basename(infile).replace(ext, "")
        m = re.match("^(.+)(" + reads[0] + "|" + reads[1] + ")$", fname)
        if m:
            bname = m.group(1)
            if bname not in infiles_dic:
                infiles_dic[bname] = [infile]
            else:
                infiles_dic[bname].append(infile)
    if not infiles_dic:
        sys.exit("Error: No fastq file has been found to be checked.")
    values_length = [len(x) for x in infiles_dic.values()]
    if min(values_length) == 2:
        pairedEnd = True
    elif min(values_length) == 1 and max(values_length) == 2:
        sys.exit("Error: The directory contains a mixture of paired-end and single-end data!")
    return pairedEnd