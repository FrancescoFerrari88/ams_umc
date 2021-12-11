#!/usr/bin/env python

# general functions for workflow #############################################
##############################################################################
import os
import yaml
import sys

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
