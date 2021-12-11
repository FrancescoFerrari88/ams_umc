#!/usr/bin/env python

# general functions for workflow #############################################
##############################################################################
import subprocess
import os
import yaml
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
